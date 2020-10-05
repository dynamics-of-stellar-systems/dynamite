
set InstallationComplete = "False"

echo ''
echo "install_mysifdec : Do you want to 'make all' in"
echo " $MYSIFDEC now (Y/n)?"
set YESNO = ""
while ( $YESNO != 'Y' && $YESNO != 'N' )
    set YESNO = $<
    if ( $YESNO == "" ) set YESNO = "Y"
    set YESNO = `echo $YESNO | tr a-z A-Z`
end

switch ( $YESNO )
    case Y: 
	$MAKE -s clean
	$MAKE -s all
	$MAKE -s clean
	set InstallationComplete = "True"
	breaksw
    case N:
	echo ' To complete the installation, type'
	echo '  make -s all'
	echo ' in the directory $MYSIFDEC'
	echo ''
	echo '  [Installation NOT complete]'
	exit 0
	breaksw
endsw

# Final environment check

if( $mysifdecIsCWD == 1 ) then
    echo "install_mysifdec : Warning : The environment variable MYSIFDEC is not set"
    echo "  It should point to your working instance of SifDec."
endif

if( ! $?MASTSIF ) then
    echo "install_mysifdec : Warning : The environment variable MASTSIF is not set"
    echo "  It should point to your main SIF repository."
endif

if( "InstallationComplete" == "True" ) then

    echo ''
    echo ' If all your environment variables are properly set,'
    echo ' you can test your installation using the command:'
    echo '   sifdecode ROSENBR'
    echo ''
    echo ' (Note that once your PATH has been updated you need'
    echo ' to use the  rehash  command to make it active.)'
    echo ''
    echo '   [Installation complete]'
    echo ''
    echo ' -------------------------------------------------------'
    echo ' '

endif

# End
