## show.awk
## Displays possible parameter values in a SIF file
##
## D. Orban, 2001, for CUTEr

BEGIN{ FS=" "; RS="\n"}{
    output = $2"="$3"\t";

    lstar = split( $1, fieldArray, "\*");
    if( lstar == 1 ){
	output = output"("$1")\t-default value-";
    }
    else if( lstar == 2 ){
	output = output"("fieldArray[2]")\t\t";
    }

    # See if some comments follow $-PARAMETER.
    # If so, print them.

    if( NF > 4 ){
	output = output"\tcomment: ";
	for ( f=5; f<=NF; f++ ){
	    output = output" "$f;
	}
    }
    else{
	output = output"\tuncommented ";
    }

    print output;
}
