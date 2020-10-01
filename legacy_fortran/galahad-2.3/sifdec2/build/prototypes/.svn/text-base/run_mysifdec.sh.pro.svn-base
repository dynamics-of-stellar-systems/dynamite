
InstallationComplete="False"

echo ''
echo "install_mysifdec : Do you want to 'make all' in"
echo " $MYSIFDEC now (Y/n)?"
YESNO=""
while [[  $YESNO != 'Y' && $YESNO != 'N'  ]]; do
    read YESNO
    [[  $YESNO == ""  ]] && YESNO="Y"
    YESNO=`echo $YESNO | tr a-z A-Z`
done

case  $YESNO  in
    Y) 
	$MAKE -s clean
	$MAKE -s all
	$MAKE -s clean
	InstallationComplete="True"
	;;
    N)
	echo ' To complete the installation, type'
	echo '  make -s all'
	echo ' in the directory $MYSIFDEC'
	echo ''
	echo '  [Installation NOT complete]'
	exit 0
	;;
esac

# Final environment check

if [[  $mysifdecIsCWD == 1  ]]; then
    echo "  Warning : The environment variable MYSIFDEC is not set"
    echo "  It should point to your working instance of SifDec."
fi

if [[ ${MASTSIF+set} != 'set'  ]]; then
    echo "  Warning : The environment variable MASTSIF is not set"
    echo "  It should point to your main SIF repository."
fi

if [[  "InstallationComplete" == "True"  ]]; then

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

fi

# End
