
# Environment check

if( ! $?SIFDEC ) then
    echo ' SIFDEC is not set.'
    echo ' Set it to the appropriate value and re-run.'
    echo ' Aborting.'
    exit 1
endif

@ mysifdecIsCWD = 0

if( -e ./install_mysifdec ) then
    if( ! $?MYSIFDEC ) then
	@ mysifdecIsCWD = 1
	setenv MYSIFDEC `dirs -l`
	setenv MYSIFDEC `echo $MYSIFDEC | $SED 's"/tmp_mnt""'`
    endif
else
    echo ' Launch install_mysifdec from its home directory.'
    echo ' Aborting.'
    exit 2
endif

# End of environment check

# Calling-sequence for install_mysifdec

set print_instructions = "false"

@ args_required = 0

if( -d $MYSIFDEC/single & -d $MYSIFDEC/double ) @ args_required = 1

if( $#argv != $args_required ) then
    set print_instructions = "true"
else if( $args_required == 1 ) then
    if( "$1" != "-DSinglePrecision" & "$1" != "-DDoublePrecision" ) then
      set print_instructions = "true"
    endif
endif

if( $print_instructions == "true" ) then
  if ( $args_required == 1 ) then
    echo "Use: install_mysifdec arg"
    echo "where arg is either"
    echo "   -DSinglePrecision  to install the single precision decoder"
    echo "   -DDoublePrecision  to install the double precision decoder"
    exit 0
  else if ( $args_required == 0 ) then
    echo " Your install_mysifdec expects no argument"
    exit 0
  endif
endif

# End of arg check

