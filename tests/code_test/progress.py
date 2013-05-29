#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

import sys

class Progress:

    def __init__(self,name,maxnum,barlength=80):
        print "%20s : " % name,
        self.maxnum = maxnum
        self.barlength = barlength
        self.count = 0
        self.steplength = self.maxnum//self.barlength

    def done(self):
        print

    def __call__(self):
        self.count += 1
        if self.count==self.steplength:
            self.count = 0
            sys.stdout.write(".")
            sys.stdout.flush()


def progress(name,seq):
    l = len(seq)
    if l==0:
        return
    p = Progress(name,l)
    for i in seq:
        p()
        yield i
    p.done()
    
