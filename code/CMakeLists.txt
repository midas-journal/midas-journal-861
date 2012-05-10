PROJECT( Mattes.Register.Series)

FIND_PACKAGE ( ITK )
IF (  ITK_FOUND  )
         INCLUDE(  ${USE_ITK_FILE}  )
ENDIF(   ITK_FOUND   )


INCLUDE_DIRECTORIES(${Mattes.Register.Series_SOURCE_DIR})

ADD_EXECUTABLE( Mattes.Register.Series  Mattes.Register.Series.cxx )

TARGET_LINK_LIBRARIES ( 
  Mattes.Register.Series 
  ITKBasicFilters 
  ITKCommon 
  ITKIO
  ITKAlgorithms
  ITKNumerics
)