# ####################################################################
# This script launch the compilation of both reports (one per sample)
# and rename them accordingly with the sample name
# ####################################################################

library( knitr)
library( rmarkdown)

### Define working folder (contains R/Rmd file for current sample, parent contains global project files)
library( funr)
if( exists( "snakemake")){
  cat("\nWorking in Snakemake mode")
  WORKING_DIR = snakemake@scriptdir
}else{
  cat("\nWorking in local script mode")
  WORKING_DIR = tryCatch(
    {
      dirname( sys.script())
    },
    error=function(cond) {
      cat("\nWorking in local R session mode")
      return( getwd())
    }
  )
}

### Load parameters
# Define an environment that will contain parameters
paramsEnv = new.env();

# Assign the WORKING_DIR to the paramsEnv
assign( "WORKING_DIR" , WORKING_DIR , env = paramsEnv )

# Load file defining global parameters
globalParamsFilePath = file.path( WORKING_DIR, "../globalParams.R");
if(file.exists(globalParamsFilePath)) {
  source( globalParamsFilePath, local = paramsEnv);
} else {
  warning("The file 'globalParamsFilePath.R' containing global parameters is missing.");
  
}

# Load file defining sample parameters
sampleParamsFilePath = file.path( WORKING_DIR, "../sampleParams.R");
if(file.exists(sampleParamsFilePath)) {
  source( sampleParamsFilePath, local = paramsEnv);
} else {
  warning("The file 'sampleParamsFilePath.R' containing sample parameters is missing.");
}

# Load file defining analysis parameters
analysisParamsFilePath = file.path( WORKING_DIR, "analysisParams.R");
if(file.exists(analysisParamsFilePath)) {
  source( analysisParamsFilePath, local = paramsEnv);
} else {
  warning("The file 'analysisParamsFilePath.R' containing analysis-specific parameters is missing.");
}

# Compile the HTML report
assign( "ORIGINAL_PATH_ANALYSIS_OUTPUT" , paramsEnv$PATH_ANALYSIS_OUTPUT , env = paramsEnv )

# for( FINDCLUSTERS_RESOLUTION in c( 0.6, 0.8, 1.0, 1.2, 1.4)){

  FINDCLUSTERS_RESOLUTION = 0.8

  assign( "FINDCLUSTERS_RESOLUTION" , FINDCLUSTERS_RESOLUTION , env = paramsEnv )
  assign( "PATH_ANALYSIS_OUTPUT" , file.path( paramsEnv$ORIGINAL_PATH_ANALYSIS_OUTPUT, paste0("Resolution_", FINDCLUSTERS_RESOLUTION)) , env = paramsEnv )
  
    # Clean the global Environment (but not the paramsEnv)
  list_variables = ls()
  list_variables = list_variables[ list_variables != "paramsEnv"]
  rm( list = list_variables)
  
  # Assign loaded values in paramsEnv to current environment
  invisible( lapply( ls( paramsEnv, all.names = TRUE), function(x, envir)
  { 
    assign( x = x, 
            value = get( x, pos = paramsEnv), 
            pos = envir)
  }, 
  environment()));
  
  # rmarkdown::render( input = file.path( WORKING_DIR, "Report.Rmd"),
  #                    output_dir = PATH_ANALYSIS_OUTPUT,
  #                    output_file  = paste0( SCIENTIFIC_PROJECT_NAME, "_", EXPERIMENT_NAME, "_", ANALYSIS_STEP_NAME, "_Sample", SAMPLE_NAME, ".html"),
  #                    quiet = FALSE)
  

# }


