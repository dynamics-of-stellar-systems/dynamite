#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

class MyRand(object):
    '''
    Class that provides the function random() which returns a 'random number'
    in the open(!) interval (0,1). The class has been created for TESTING and
    DEVELOPMENT purposes and ALWAYS PRODUCES THE SAME SEQUENCE of 'random
    numbers' (as long as no multiple threads are running). This also works
    'cross-language' with the Fortran implementation. If desired, the sequence
    can be changed by assigning the class attribute MyRand.idum a different
    negative integer.
    '''

    idum = -42
    iv = None
    iy = None

    def __init__(self):
        self.idum = self.__class__.idum
        self.iv = self.__class__.iv
        self.iy = self.__class__.iy

    def random(self):
        '''
        From Numerical Recipes in F77, 2nd. Edition, corresponds to ran1.
        “Minimal” random number generator of Park and Miller with Bays-Durham
        shuffle and added safeguards. Returns a uniform random deviate between
        0.0 and 1.0 (exclusive of the endpoint values). Call with self.idum a
        negative integer to initialize; thereafter, do not alter idum between
        successive deviates in a sequence. RNMX should approximate the largest
        floating value that is less than 1.

        Returns
        -------
        np.float32
            Next 'random number' in the sequence.

        '''
        # INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
        # REAL ran1,AM,EPS,RNMX
        # IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836
        # NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS
        IA=16807
        IM=2147483647
        AM=np.float32(1./IM)
        IQ=127773
        IR=2836
        NTAB=32
        NDIV=1+int((IM-1)/NTAB)
        EPS=np.float32(1.2e-7)
        RNMX=np.float32(1.-EPS)
    
        # INTEGER j,k,iv(NTAB),iy
        # SAVE iv,iy
        # DATA iv /NTAB*0/, iy /0/
        if self.iv is None:
            self.iv = np.int_([0]*NTAB)
        if self.iy is None:
            self.iy = 0
        
        if self.idum <= 0 or self.iy == 0:
            self.idum=max(-self.idum,1)
            # do 11 j=NTAB+8,1,-1
            for j in range(NTAB+8,0,-1):
                # print(f'\t1 ran1.idum={ran1.idum}')
                k = self.idum // IQ
                self.idum = IA*(self.idum-k*IQ)-IR*k
                if self.idum < 0:
                    self.idum += IM
                if j <= NTAB:
                    self.iv[j-1] = self.idum
            self.iy = self.iv[0]
        # print(f'\t2 ran1.idum={ran1.idum}')
        k = self.idum // IQ
        self.idum = IA*(self.idum-k*IQ)-IR*k
        if self.idum < 0:
            self.idum += IM
        # print(f'\t3 ran1.iy={ran1.iy}')
        j = 1 + self.iy // NDIV
        self.iy = self.iv[j-1]
        self.iv[j-1] = self.idum
        return min(np.float32(AM*self.iy),RNMX)
