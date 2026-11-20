/*==============================================================================
  Program   : local_config_template.sas
  Project   : synpuf-survival-sas
  Purpose   : Template for the one setting that is specific to your computer.
  Usage     : Copy this file to sas/local_config.sas in the same folder and
              change the path below to the folder that holds your unzipped
              CMS DE-SynPUF Sample 1 Beneficiary Summary Files.
              00_config.sas reads sas/local_config.sas automatically.
              .gitignore excludes sas/local_config.sas, so your path is never
              committed. As an alternative, set the SYNPUF_DATA_ROOT
              environment variable and leave local_config.sas absent.
==============================================================================*/
%let data_root = C:\path\to\sample_01;
