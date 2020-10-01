
# Environment check

if( ! $?CUTER ) then
    echo ' CUTER is not set.'
    echo ' Set it to the appropriate value and re-run.'
    echo ' Aborting.'
    exit 1
endif

@ mycuterIsCWD = 0

if( -e ./install_mycuter ) then
    setenv MYCUTER `dirs -l`
    setenv MYCUTER `echo $MYCUTER | $SED 's"/tmp_mnt""'`
else
    echo ' Launch install_mycuter from its home directory.'
    echo ' Aborting.'
    exit 2
endif

# End of environment check

# Calling sequence for install_mycuter

set print_instructions = "false"

@ args_required = 0

if( -d $MYCUTER/single & -d $MYCUTER/double ) @ args_required = 1

if( $#argv != $args_required ) then
    set print_instructions = "true"
else if( $args_required == 1 ) then
    if( "$1" != "-DSinglePrecision" & "$1" != "-DDoublePrecision" ) then
      set print_instructions = "true"
    endif
endif

if( $print_instructions == "true" ) then
  if ( $args_required == 1 ) then
    echo "Use: install_mycuter arg"
    echo "where arg is either"
    echo "   -DSinglePrecision  to install the single precision tools"
    echo "   -DDoublePrecision  to install the double precision tools"
    exit 0
  else if ( $args_required == 0 ) then
    echo " Your install_mycuter expects no argument"
    exit 0
  endif
endif

# End of arg check
