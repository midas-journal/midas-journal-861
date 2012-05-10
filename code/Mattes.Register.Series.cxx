/*=========================================================================

Author: Ken Urish
Program: Register Large Serie
Date: 20 April 2011

Objective:
based on Multimodality Mattes,
Uses the crop image to select region to register
Input is a text file with a list of file names to register
NOTE BUG ALERT
Make sure you make sure the order of DESS and T2 file names match the order that they are read in.


/***********************************************************************************


Base
100822_Register_MattesRegister
VerserRigid3DTransform project
ImageRegistrationExample7.cxx - Mattes Metric
ImageRegistrationExample8.cxx - Verser Transform

Input:
-Directory: "Fixed" of dicom images
-Directory: "Moving" of dicom images

Output:
-Txt file of iteration results:
-Output.mhd
-DifferenceBefore.mhd
-DifferenceAfter.mhd

To do: 
-Have these sent out as dicoms. Thus can easily view slice by slice.

*************************************************************************************/

//Theory Note
//  The parameter space of the \code{VersorRigid3DTransform} is not a vector
//  space, due to the fact that addition is not a closed operation in the space
//  of versor components. This precludes the use of standard gradient descent
//  algorithms for optimizing the parameter space of this transform. A special
//  optimizer should be used in this registration configuration. The optimizer
//  designed for this transform is the VersorRigid3DTransformOptimizer}. This optimizer uses Versor
//  composition for updating the first three components of the parameters
//  array, and Vector addition for updating the last three components of the
//  parameters array~\cite{Hamilton1866,Joly1905}.


//  The metric requires two parameters: number of entropy bins (50 in 2d) 
//  number of spatial samples to compute the density estimates. 
//  Number of bins may have dramatic effects on the optimizer's behavior. 
//  Number of spatial samples to be used depends on the content of the image. 
//  If the images are smooth and do not contain much detail, then using approximately $1$ percent of the
//  pixels will do. On the other hand, if the images are detailed, it may be
//  necessary to use a much higher proportion, such as $20$ percent.
  
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


//ITK Registration
#include "itkImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkCheckerBoardImageFilter.h"
#include "itkFixedArray.h"

//ITK Filters
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
        
//ITK Basic
#include "itkCommand.h" //for command oberver
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"



//Steak and Potatoes
//Im not using all of these
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <itksys/SystemTools.hxx>
#include <time.h> //to time program


//----------------------------------------------------------------------------------
// COMAND OBSERVER CLASS
//----------------------------------------------------------------------------------
class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  typedef itk::VersorRigid3DTransformOptimizer     OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer = 
        dynamic_cast< OptimizerPointer >( object );
      if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
        return;
        }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue() << "   ";
      std::cout << optimizer->GetCurrentPosition() << std::endl;

	  /*
	  std::fstream outputFile; 
      outputFile.open("output.txt",std::ios::out | std::ios::app);
	  outputFile <<optimizer->GetCurrentIteration()<<"\t";
	  outputFile <<optimizer->GetValue()<<"\t";
	  outputFile <<optimizer->GetCurrentPosition()<<std::endl;
      outputFile.close();
	  */
    }
};


//-----------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// ITK to VTK
// Connect the given itk::VTKImageExport filter to given vtkImageImport filter.
//-------------------------------------------------------------------------------
template <typename ITK_Exporter, typename VTK_Importer>
void ConnectPipelines(ITK_Exporter exporter, VTK_Importer* importer)
{
  importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
  importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
  importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
  importer->SetSpacingCallback(exporter->GetSpacingCallback());
  importer->SetOriginCallback(exporter->GetOriginCallback());
  importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
  importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
  importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
  importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
  importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
  importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
  importer->SetCallbackUserData(exporter->GetCallbackUserData());
}






int main( int argc, char *argv[] )
{

//----------------------------------------------------------------------------------
// INPUT PARAMETERS
//----------------------------------------------------------------------------------

  // Read input file
  std::string inputFileName = argv[1]; 
  std::fstream inputFile;
  inputFile.open(inputFileName.c_str(),std::ios::in);
  inputFile.seekg(0, std::ios::beg);	
  inputFile.clear();

  //Optimizer Parameters
  const double maxStepLength = atof (argv[2]); //0.2000; 
  const double minStepLength = atof (argv[3]); //0.0001;
  const unsigned int iterations = atoi (argv[4]); //10;
  
  //Mattes Metric Parameters
  const double metricPercent = atof (argv[5]); //0.01 to 0.2
  unsigned int numberOfBins = atoi (argv[6]); //50;
  const double translationScaleFactor = atof (argv[7]); //10.0;
 
  //This figures out which parts of image to crop to get even overlap
  double regionPercent = atof (argv[8]); //0.05;
  int numSlices = atoi (argv[9]); //2;

  //checker board layout output
  unsigned int xCheckerNum = atoi (argv[10]);
  unsigned int yCheckerNum = atoi (argv[11]);
  unsigned int zCheckerNum = atoi (argv[12]);
  
  //
  unsigned int fixedMinValue = atoi (argv[13]);
  unsigned int fixedMaxValue = atoi (argv[14]);
  unsigned int movingMinValue = atoi (argv[15]);
  unsigned int movingMaxValue = atoi (argv[16]);

  //Start of big loop for rest of program to read filenames
  std::string fixedDirectoryName;
  std::string movingDirectoryName;
  


  //Opening Observer File
    const char * txtFile = "output.111228.txt";
  std::fstream outputFile; 
  outputFile.open(txtFile,std::ios::out);
  outputFile <<"maxStepLength (norm 0.2): "<<maxStepLength<<std::endl;
  outputFile <<"minStepLength (norm 0.0001): "<<minStepLength<<std::endl;
  outputFile <<"iterations (norm 200): "<<iterations<<std::endl;
  outputFile <<"metric percent (norm 0.25): "<<metricPercent<<std::endl;
  outputFile <<"number of bins (norm 150): "<<numberOfBins<<std::endl;
  outputFile <<"translation scale factor (norm 10): "<<translationScaleFactor<<std::endl;
  outputFile <<"region percent (norm 0-0.5): "<<regionPercent<<std::endl;
  outputFile <<"number of z-slices (norm 0-2): "<<numSlices<<std::endl<<std::endl<<std::endl<< std::endl<< std::endl;
  outputFile <<"id:\t  value:\t versorX:\t versorY:\t magnitudeRotation:\t versorZ:\t changeX:\t changeY:\t changeZ:\t magnitudeTranslation:\t iterations:\t metric:\t time(s):"<< std::endl;  
  outputFile.close();

  while (!inputFile.eof())
  {

  //----------------------------------------------------------------------------------
  // TIMER
  //----------------------------------------------------------------------------------  
  time_t startTime;
  time_t stopTime;
  startTime = time (NULL);

  //In OAI, DESS is listed first; DESS is the moving as has higher resolution
  inputFile >> movingDirectoryName >> fixedDirectoryName;

  std::string fixedDirectory  = fixedDirectoryName;  //T2
  std::string movingDirectory = movingDirectoryName;  //DESS
  
  
  std::string outputDirectory = fixedDirectory + "/" + "output_111228";  //Transformed moving image (T2)
  std::string afterDirectory  = fixedDirectory + "/" + "after_111228";  //difference after
  std::string textFile        = fixedDirectory + "/" + "output_111228.txt";  //difference after
  
  itksys::SystemTools::MakeDirectory(afterDirectory.c_str());
  itksys::SystemTools::MakeDirectory(outputDirectory.c_str());
  
  std::cout<< "Fixed: " << fixedDirectory.c_str()<<std::endl;
  std::cout<< "Moving: " << movingDirectory.c_str()<<std::endl;
  std::cout<< afterDirectory.c_str()<<std::endl;
  std::cout<< textFile.c_str()<<std::endl;
  


//----------------------------------------------------------------------------------
// Keg-o-Typedefs
//----------------------------------------------------------------------------------
  const    unsigned int    Dimension = 3;
  typedef  signed short   InputPixelType;
  typedef  signed short   OutputPixelType;
  typedef  unsigned char   vtkPixelType;

  typedef itk::Image< InputPixelType, Dimension >  FixedImageType;
  typedef itk::Image< InputPixelType, Dimension >  MovingImageType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::Image< OutputPixelType, 2 > OutputSliceType;

  typedef itk::Image< vtkPixelType, Dimension > vtkImageType;


//----------------------------------------------------------------------------------
// DICOM SERIES READER
//----------------------------------------------------------------------------------
  typedef itk::ImageSeriesReader< FixedImageType > FixedReaderType;
  typedef itk::ImageSeriesReader< MovingImageType > MovingReaderType;
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  typedef itk::GDCMImageIO  ImageIOType;
  
  ImageIOType::Pointer gdcmIO = ImageIOType::New();
  FixedReaderType::Pointer fixedReader  = FixedReaderType::New();
  MovingReaderType::Pointer movingReader  = MovingReaderType::New();
  NamesGeneratorType::Pointer fixedNamesGenerator = NamesGeneratorType::New();
  NamesGeneratorType::Pointer movingNamesGenerator = NamesGeneratorType::New();
  
  //reading fixed image files
  fixedNamesGenerator->SetInputDirectory( fixedDirectory );
  const FixedReaderType::FileNamesContainer & fixedFilenames = fixedNamesGenerator->GetInputFileNames();
  fixedReader->SetImageIO( gdcmIO );
  

  //This section will read in T2 as new image. Will uso dicom tags to screen only first echo.
  //Only pass the first echo of T2 file names (ie 27; not all 196)  
  //(0018,0081) EchoTime needs to be around 10ms


  FixedReaderType::Pointer sortFixedReader  = FixedReaderType::New();
  NamesGeneratorType::Pointer sortFixedNamesGenerator = NamesGeneratorType::New();
  sortFixedNamesGenerator->SetUseSeriesDetails( true );
  sortFixedNamesGenerator->AddSeriesRestriction("0018|0081" );
  sortFixedNamesGenerator->SetDirectory( fixedDirectory );

  typedef std::vector< std::string >    SeriesIdContainer;
    
  const SeriesIdContainer & seriesUID = sortFixedNamesGenerator->GetSeriesUIDs();  
  SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
  SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

  while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      seriesItr++;
      }

  //this uses the first series identifier as a default for selecting the series we will use for T2 images
  //an easy upgrade is adding argument where user can define their own series.
  std::string seriesIdentifier;
  seriesIdentifier = seriesUID.begin()->c_str();

  std::cout << std::endl << std::endl;
  std::cout << "Now reading series: " << std::endl << std::endl;
  std::cout << seriesIdentifier << std::endl;

  typedef std::vector< std::string >   FileNamesContainer;
  FileNamesContainer fixedFileNames;
  fixedFileNames = sortFixedNamesGenerator->GetFileNames( seriesIdentifier );
  
  unsigned int numSortFixedFilenames =  fixedFileNames.size();
  std::cout << numSortFixedFilenames << std::endl;
  for(unsigned int fni = 0; fni<numSortFixedFilenames; fni++)
  { std::cout << "fixed file names: "<< fixedFileNames[fni] << std::endl; }  
  
  fixedReader->SetFileNames( fixedFileNames );


  //reading moving image files
  movingNamesGenerator->SetInputDirectory( movingDirectory );
  const MovingReaderType::FileNamesContainer & movingFilenames = movingNamesGenerator->GetInputFileNames();
  movingReader->SetImageIO( gdcmIO );
  movingReader->SetFileNames( movingFilenames );
  unsigned int numMovingFilenames =  movingFilenames.size();
  std::cout << "moving file names number: "<< numMovingFilenames << std::endl; 
  
  //need update so register method can get region size on call
  try { fixedReader->Update(); } 
  catch( itk::ExceptionObject & err ) 
    { 
	std::cerr << "ITK says: What? You cant even read an image!" << std::endl; 
    std::cerr << err << std::endl;
	getch();
    return EXIT_FAILURE;
    } 

 //need update so register method can get region size on call
  try { movingReader->Update(); } 
  catch( itk::ExceptionObject & err ) 
    { 
	std::cerr << "ITK says: What? You cant even read an image!" << std::endl; 
    std::cerr << err << std::endl;
	getch();
    return EXIT_FAILURE;
    } 


//----------------------------------------------------------------------------------
// Declarations for Registration
//----------------------------------------------------------------------------------
  
  //Putting this before vtk gets rid of some bugs

  FixedImageType::Pointer fixedImagePtr = fixedReader->GetOutput();
  FixedImageType::SpacingType spacing = fixedImagePtr->GetSpacing();
  FixedImageType::SizeType    size    = fixedImagePtr->GetLargestPossibleRegion().GetSize();
  FixedImageType::PointType   origin  = fixedImagePtr->GetOrigin();

  MovingImageType::Pointer movingImagePtr;

  typedef itk::VersorRigid3DTransform< double > TransformType;
  typedef itk::VersorRigid3DTransformOptimizer  OptimizerType;
  typedef itk::LinearInterpolateImageFunction< MovingImageType, double> InterpolatorType;
  typedef itk::ImageRegistrationMethod<FixedImageType, MovingImageType>  RegistrationType;
  typedef itk::MattesMutualInformationImageToImageMetric< FixedImageType, MovingImageType > MetricType;

  //The components of the registration process
  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
  TransformType::Pointer  transform = TransformType::New();

  fixedImagePtr = fixedReader->GetOutput();
  movingImagePtr = movingReader->GetOutput();
  
  FixedImageType::RegionType regionRegister;
  FixedImageType::SizeType sizeRegister;
  FixedImageType::IndexType indexRegister;



//----------------------------------------------------------------------------------
// VERSER REGISTRATION
//----------------------------------------------------------------------------------

  // SELECT REGION
  // Sets up the index and the size of the region to register

  // Begeining of '1'. Gross region of interest. 
  //Index of Registration Region (IN PIXELS)
  indexRegister[0] = regionPercent*size[0]; 
  indexRegister[1] = regionPercent*size[1];
  indexRegister[2] = numSlices; //z can be less resolution

  //Size of registration Region (IN PIXELS)
  sizeRegister[0] = size[0] - regionPercent*size[0]; 
  sizeRegister[1] = size[1] - regionPercent*size[1];
  sizeRegister[2] = size[2] - numSlices;

  //end of gross region of interest

  regionRegister.SetIndex(indexRegister);
  regionRegister.SetSize (sizeRegister);
  
  //Check region dimensions
/*
  outputFile.open(txtFile,std::ios::out | std::ios::app);
  outputFile <<"Region registered: "<<regionRegister<<std::endl;
  std::cout <<"Region registered: "<<regionRegister<<std::endl;
  outputFile.close();
*/

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetInterpolator(  interpolator  );    
  registration->SetFixedImage(    fixedReader->GetOutput()    );
  registration->SetMovingImage(   movingReader->GetOutput()   );
  registration->SetTransform(     transform );
  registration->SetFixedImageRegion(regionRegister);

  // Metric Type
  // One mechanism for bringing the Metric to its limit is to disable the
  // sampling and use all the pixels present in the FixedImageRegion. 
  unsigned int numberOfSamples = 10000;
  const unsigned long numberOfImagePixels = 
					fixedImagePtr->GetLargestPossibleRegion().GetNumberOfPixels();
  const unsigned long numberOfSpatialSamples =
					static_cast<unsigned long>( numberOfImagePixels * metricPercent );
  metric->SetNumberOfSpatialSamples( numberOfSpatialSamples );
  metric->SetNumberOfHistogramBins( numberOfBins );


  //Transform Initializer
  //Not necessary to call Update() on reader since the 
  //CenteredTransformInitializer will do it as part of its computations. 
  //If want geometrical centers call: GeometryOn()} 
  //If want mass centers call: MomentsOn().
  //Computation of the center and translation is triggered by InitializeTransform()
  typedef itk::CenteredTransformInitializer< TransformType, FixedImageType, 
	                                         MovingImageType>  TransformInitializerType;
  TransformInitializerType::Pointer initializer = TransformInitializerType::New();
  initializer->SetTransform(   transform );
  initializer->SetFixedImage(  fixedReader->GetOutput() );
  initializer->SetMovingImage( movingReader->GetOutput() );
  initializer->MomentsOn();
  initializer->InitializeTransform();

  //  The rotation part of the transform is initialized using a versor, which is simply a unit quaternion.  
  //  The versor itself defines the type of the vector used to indicate the rotation axis.
  //  This creates a versor object and initialize its parameters by passing a
  //  rotation axis and an angle.

  typedef TransformType::VersorType  VersorType;
  typedef VersorType::VectorType     VectorType;
  
  VersorType     rotation;
  VectorType     axis;
  axis[0] = 0.0;
  axis[1] = 0.0;
  axis[2] = 1.0;

  const double angle = 0;
  rotation.Set(  axis, angle  );
  transform->SetRotation( rotation );
  
  //  We now pass the parameters of the current transform as the initial
  //  parameters to be used when the registration process starts. 
  registration->SetInitialTransformParameters( transform->GetParameters() );
  
  //  Another significant difference in the metric is that it computes the
  //  negative mutual information and hence we need to minimize the cost
  //  function in this case.
  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
  
  optimizerScales[0] = 1.0;
  optimizerScales[1] = 1.0;
  optimizerScales[2] = 1.0;
  optimizerScales[3] = 1.0 / ( translationScaleFactor*spacing[0]*size[0] );
  optimizerScales[4] = 1.0 / ( translationScaleFactor*spacing[1]*size[1] );
  optimizerScales[5] = 1.0 / ( translationScaleFactor*spacing[2]*size[2] );
  
  optimizer->SetScales( optimizerScales );
  optimizer->SetMaximumStepLength(  maxStepLength);  
  optimizer->SetMinimumStepLength(  minStepLength); 
  optimizer->SetNumberOfIterations( iterations);
  optimizer->MaximizeOff();

  // Create the Command observer and register it with the optimizer.
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer);

  //Lets see if we can register
  try { registration->StartRegistration(); } 
  catch( itk::ExceptionObject & err ) 
    { 
	std::cerr << "ITK says: You wrong - you will respect my authority!" << std::endl; 
    std::cerr << err << std::endl; 
	getch();
    return EXIT_FAILURE;
    } 

//----------------------------------------------------------------------------------
// TIMER
//----------------------------------------------------------------------------------
  stopTime = time (NULL);
  stopTime = (stopTime-startTime); 


//----------------------------------------------------------------------------------
// NUMERICAL OUTPUT
//----------------------------------------------------------------------------------
  OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
  const double versorX              = finalParameters[0];
  const double versorY              = finalParameters[1];
  const double versorZ              = finalParameters[2];
  const double finalTranslationX    = finalParameters[3];
  const double finalTranslationY    = finalParameters[4];
  const double finalTranslationZ    = finalParameters[5];

  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  const double bestValue = optimizer->GetValue();

  transform->SetParameters( finalParameters );
  TransformType::MatrixType matrix = transform->GetRotationMatrix();
  TransformType::OffsetType offset = transform->GetOffset();


  // Print out results
  std::cout << std::endl << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " versor X      = " << versorX  << std::endl;
  std::cout << " versor Y      = " << versorY  << std::endl;
  std::cout << " versor Z      = " << versorZ  << std::endl;
  std::cout << " Translation X = " << finalTranslationX  << std::endl;
  std::cout << " Translation Y = " << finalTranslationY  << std::endl;
  std::cout << " Translation Z = " << finalTranslationZ  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Best Metric   = " << bestValue          << std::endl<< std::endl;
  std::cout << " Matrix = " << std::endl << matrix << std::endl;
  std::cout << " Offset = " << std::endl << offset << std::endl;

  //The versor can be converted to the magnitude of the rotation.
  //See Louis email post: http://www.itk.org/pipermail/insight-users/2009-March/029584.html
  //The magntiude is in radians and needs ocnverted to degrees.
  double radiansToDegrees = 45.0/atan(1.0);
  double magnitudeRotation =  2.0 * asin(sqrt(versorX*versorX+ versorY*versorY+ versorZ*versorZ))*radiansToDegrees;

  //magnitude of the final translation in pixels
  double magnitudeTranslation =  sqrt(finalTranslationX/spacing[0]*finalTranslationX/spacing[0]+finalTranslationY/spacing[1]*finalTranslationY/spacing[1]+finalTranslationZ/spacing[2]+finalTranslationZ/spacing[2]);

  //Opening Observer File
  outputFile.open(txtFile,std::ios::out | std::ios::app);
  outputFile << fixedDirectory.c_str()<<"\t"<< versorX <<"\t"<< versorY <<"\t"<< versorZ <<"\t"<< magnitudeRotation <<"\t"<< finalTranslationX <<"\t"<< finalTranslationY <<"\t"<< finalTranslationZ <<"\t"<< magnitudeTranslation <<"\t"<< numberOfIterations <<"\t"<< bestValue <<"\t"<< stopTime<< std::endl;
  outputFile.close();

  //----------------------------------------------------------------------------------
  // TRANSFORMS AND RESAMPLING
  //----------------------------------------------------------------------------------
  TransformType::Pointer finalTransform = TransformType::New();
  finalTransform->SetCenter( transform->GetCenter() );
  finalTransform->SetParameters( finalParameters );

  typedef itk::ResampleImageFilter<MovingImageType,FixedImageType > ResampleFilterType;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  FixedImageType::Pointer fixedImage = fixedImagePtr;
  resampler->SetTransform( finalTransform );
  resampler->SetInput( movingImagePtr );
  resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin( fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  //resampler->SetOutputDirection( fixedImage->GetDirection() );
  //resampler->SetDefaultPixelValue( 100 );


 
//----------------------------------------------------------------------------------
// WRITERS
//----------------------------------------------------------------------------------
  
  //Series writer declaration for After and Before shots  
  typedef itk::CastImageFilter<FixedImageType, OutputImageType > CastFilterType;
  typedef itk::ImageSeriesWriter< OutputImageType, OutputSliceType > SeriesWriterType;

  CastFilterType::Pointer  caster =  CastFilterType::New();
  caster->SetInput( resampler->GetOutput() );

  SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
  sortFixedNamesGenerator->SetOutputDirectory(outputDirectory.c_str());    
  //seriesWriter->SetInput(caster->GetOutput() );
  seriesWriter->SetInput(resampler->GetOutput() );
  seriesWriter->SetImageIO( gdcmIO );
  seriesWriter->SetFileNames( sortFixedNamesGenerator->GetOutputFileNames() );
  
  try { seriesWriter->Update(); } 
  catch( itk::ExceptionObject & err ) 
    { 
	std::cerr << "ITK says: You wrong!" << std::endl; 
    std::cerr << err << std::endl; 
	getch();
    return EXIT_FAILURE;
    } 


	
// Generate checkerboards after registration
//----------------------------------------------------------------------------------
// ADJUST INPUT IMAGE SEQUENCE CONTRASTS
//----------------------------------------------------------------------------------

  //The following is a patch for image output. For the checker board layout, the images are easier to view
  //if there contrast/brightness is altered. The following section is used to alter the histogram and thus
  //adjust teh contrast for imrpoved visualization. 

  //Bug Note: This means the output is not the same image intensity as the original for the DESS!
  //This is being done at the end of he program so it doesnt mess with the pixel values of any of the other input images.

  //Here, we will resample the DESS to match the spectrum of the T2
/*
  typedef itk::RescaleIntensityImageFilter< FixedImageType, OutputImageType >   FixedCastFilterType;
  FixedCastFilterType::Pointer fixedCast = FixedCastFilterType::New();
  fixedCast->SetInput(fixedReader->GetOutput());
  fixedCast->SetOutputMaximum (fixedMaxValue);
  fixedCast->SetOutputMinimum (fixedMinValue);
  fixedCast->Update();  
  fixedImagePtr = fixedCast->GetOutput();
*/

  typedef itk::RescaleIntensityImageFilter< FixedImageType, OutputImageType >   MovingCastFilterType;
  MovingCastFilterType::Pointer movingCast = MovingCastFilterType::New();
//  movingCast->SetOutputMaximum (movingMaxValue);
//  movingCast->SetOutputMinimum (movingMinValue);
  
   movingCast->SetOutputMaximum (1000);
   movingCast->SetOutputMinimum (0);
  
  movingCast->SetInput(resampler->GetOutput());
  movingCast->Update(); 
  movingImagePtr = movingCast->GetOutput();


  typedef itk::FixedArray< unsigned char, 3 > ArrayType;
  ArrayType checkerNum; 
  checkerNum[0] = xCheckerNum;
  checkerNum[1] = yCheckerNum;
  checkerNum[2] = zCheckerNum;

  typedef itk::CheckerBoardImageFilter< FixedImageType > CheckerBoardFilterType;
  CheckerBoardFilterType::Pointer checker = CheckerBoardFilterType::New();
  checker->SetInput1( fixedImage );
  checker->SetInput2(  movingImagePtr );
  checker->SetCheckerPattern (checkerNum);

  //This was the old code.
  //checker->SetInput2( resampler->GetOutput() );
  //caster->SetInput( checker->GetOutput() );
  //resampler->SetDefaultPixelValue( 0 );

  // After registration
  seriesWriter->SetInput(checker->GetOutput() );
  sortFixedNamesGenerator->SetOutputDirectory(afterDirectory.c_str());
  seriesWriter->SetFileNames( sortFixedNamesGenerator->GetOutputFileNames() );
  
  try { seriesWriter->Update(); } 
  catch( itk::ExceptionObject & err ) 
    { 
	std::cerr << "ITK says: You wrong!" << std::endl; 
    std::cerr << err << std::endl; 
	getch();
    return EXIT_FAILURE;
    } 
 	
/*
  //COPY FILE TO STRING
  const char * cloneFile = "copy.txt";
  std::ifstream in(txtFile); //This is the temp file of output
  std::ofstream out(textFile.c_str());  //This is the file that is put in basde directory
  std::string s;
  while (getline(in,s))
	  out<<s<<"\n";
*/ 

} //End of While loop to read filenames

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Thanks Ken. All Done"<< std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Support medical image analysis in your local community."<< std::endl;
  std::cout << "Open source is open science"<< std::endl;

  //Well done.
  return 0;
}

