## param.awk
## Finds out whether the given parameter may be set
## to the given value within a SIF file.
##
## param.awk receives the parameter name and its value
## as arguments, in the variables pname and pval respectively.
## It also receives the 'doesmatch' variable. If doesmatch=0,
## param.awk returns all the line numbers where there is a match
## for pname but *not* for pval. Otherwise, it returns the line
## number where there is a match for both.
##
## D. Orban, 2001, for CUTEr

BEGIN{ FS=" "; RS="\n"}{

# The lines we receive have the form
# n: type parameterName     parameterValue
# or
# n:*type parameterName     parameterValue
# where n is a line number, type is IE, RE, IA,...
# parameterName is a string (which could match pname)
# and parameterValue is a value (which could match pval).
#
# In the first case, NF=4, and in the second case, NF=3.

    if( NF == 3 ){
	pnameField = 2;
	pvalField  = 3;
    }
    else{
	pnameField = 3;
	pvalField  = 4;
    }

    nonMatchingLines = "";

    if( doesmatch != 0 ){
	if( pname == $pnameField && pval == $pvalField ){
# If we have a match, print the line number and exit
	    l = split( $0, fieldArray, ":" );
	    print fieldArray[1];
	    next;
	}
    }
    else{
	if( pname == $pnameField && pval != $pvalField ){
# If we have a match, print the line number and exit
	    ll = split( $0, nfieldArray, ":" );
	    nonMatchingLines = nonMatchingLines nfieldArray[1];
	}
    }

    print nonMatchingLines;
}

