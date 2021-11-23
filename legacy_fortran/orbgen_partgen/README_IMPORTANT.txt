PLEASE NOTE:

Both orbgen.f90 and partgen.f90 are legacy files which got the orbit mirroring fix but are UNTESTED.
orbgen.f90: the fix is in subroutine project_part(), original version commented out
partgen.f90: the fix is in subroutine project_part(), original version commented out

Before using them, please double-check that the signs are correct!

If you want to build partgen, add the following recipe to Makefile:

partgen : $(ORBLIBOBJ) partgen.f90
	$(link) partgen  $(ORBLIBOBJ) partgen.f90

