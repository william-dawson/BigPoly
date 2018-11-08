!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module process input based on a yaml specification file.
MODULE AutoYamlArgsModule
  !! BigDFT Modules
  USE dictionaries
  USE futile
  USE yaml_output, ONLY : yaml_dict_dump_all
  USE yaml_parse, ONLY : yaml_cl_parse, yaml_cl_parse_null, &
       & yaml_parse_from_file, yaml_parse_from_string
  USE DatabaseModule, ONLY : datastream
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: GetArgs
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE GetArgs(options, iproc)
    !> The dictionary specifying the parameters of the calculation.
    TYPE(dictionary), POINTER, INTENT(OUT) :: options
    !> The process on which this function is being run.
    INTEGER, INTENT(IN) :: Iproc
    !! Local data
    TYPE(yaml_cl_parse) :: parser
    TYPE(dictionary), POINTER :: spec
    !! Temporary Variables
    TYPE(dictionary), POINTER :: dtmp, it

    !! Setup the initial spec from file.
    CALL yaml_parse_from_string(dtmp, datastream)
    spec => dtmp .pop. 0
    CALL dict_free(dtmp)

    !! Parse the input.
    parser=yaml_cl_parse_null()
    it => dict_iter(spec)
    DO WHILE(ASSOCIATED(it))
       CALL yaml_cl_parse_option(parser, TRIM(dict_key(it)), &
            & TRIM(dict_value(it//'default')), &
            & TRIM(dict_value(it//'shorthelp')))
       it => dict_next(it)
    END DO
    CALL yaml_cl_parse_cmd_line(parser,args=options)
    CALL yaml_cl_parse_free(parser)
    IF (iproc .EQ. 0) THEN
       CALL yaml_map('Parsed options',options)
    END IF
  END SUBROUTINE GetArgs
END MODULE AutoYamlArgsModule
