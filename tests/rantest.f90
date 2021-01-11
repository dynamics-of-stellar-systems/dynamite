program rantest
implicit none
integer idum, i
double precision ran1
idum = -4242
print *, ran1(idum)
do i=1,9
  print *, ran1(1)
end do
end program rantest
