#!/usr/bin/python
'''

buildCstruct 1.2

Authors: Andrea Minoia, Martin Voegele
Date: January 27th 2018
********************************************
Description:
    Build allotropes structures of carbon, namely graphite sheets and
    carbon nanotubes. Non-periodic structure are saturated with hydrogens.
    Functional groups (OH, COOH, COO-) can be added to CNTs. 


Syntax:

    buildCstruct.py [options] outfile

    The available options are:
    --version

        shows program's version number and exit.

    -c, --credits

        display credits.

    -s, --structure

        specify the kind of structure to build. Valid structures are: armcnt, zigzagcnt and hopg.

    -p, --periodicity

        build a periodic structure for TINKER.

    -g, --geometry

        specify the geometry of the structure. For CNTs, use -g index_n cnt_length while
        for hopg use -g size_x size_y (see examples below).
        Size_x, size_y and cnt_length are in Angstom.

    -f, --functionalization

        specify the functionalization (only for CNT since version 1.2)

    --xyz

        save structure in XYZ format (version >= 1.1 only).

    --gro

        save structure in gromacs GRO format (version >= 1.1 only).

    --mol2

        save structure in mol2 format (version >= 1.2 only).

    outfile
    is the name of the file where to save the structure.


Requirements:
    Require numpy installed


Known issues and limitations:
    1) Build single wall CNT and rectangular slab of graphtie HOPG only.
    2) Build armchair and zigzag CNT only


Changelog:

    + v 1.2 - January 2018:
    	- Modified H adding, so there are no clashes between hydrogens.
    	- Added support for mol2 format
        - Added partial charges
        - Added COOH, COO- and OH groups for CNTs (ooption -f)

    + v 1.1 - March 2010:
        - improved connectivity search algorithm 
        - added support for XYZ format
        - added support for gromacs GRO format
        
    + v 1.0 - April 2009:
        - first release of buildCstruct




License:
    Freeware. You can use, modify and redistribute the source.
      
    buildCstruct is part of Yasc, a collection of freeware software, therefore it
    cames with _ABSOLUTE_ _NO_(as in NOTHING, NADA, NICHTS, NIENTE, RIEN) WARRANTY.
    The authors of the scripts collected in Yasc are _NOT_ responsible if these software
    will erase your hard disks, empty your bank account, stole your car, seduce
    your wife, shave your dog or make any kind of mess and damages, including loss
    of data or worst. By using Yasc, you _ACCEPT_ these terms.
    Once again: use at your own risk!


Contacts:
    You liked this program? You have suggestions?
    Drop me a mail at minoiaa_at_gmail.com
    Don't forget to visit the wiki: http://chembytes.wikidot.com
'''


#import modules
from sys import argv, exit,path, version

from optparse import OptionParser as OP
try:
    from numpy import zeros, pi, sin, cos, modf, ceil, sqrt
except:
    print "Numpy not installed or not in python path. I give up..."
    exit(10)
from string import lower
from os import path as path_os
from os import system


'''
#===============================================================================
#                               SUBROUTINES
#===============================================================================
'''

def getdist(at1,at2):
    ''' Calculate distance between two particles
    '''
    dist_at=sqrt((at2[0]-at1[0])**2+(at2[1]-at1[1])**2+(at2[2]-at1[2])**2)
    return dist_at

def filecheck(file):
    ''' Check if infile exists
    '''
    if path_os.isfile(file) == 0:
        found=False
    else:
        found=True
    return found

def backup_file(file):
    '''check if file exists. if not open file
    for writing, otherwhise backup the old one in
    #infile_x# with x progressive number
    '''
    tmpvar=file
    count=0
    while 1:
        found=filecheck(file)
        if found:
            count+=1
            file=tmpvar+'_bak-'+str(count)
        else:
            break
    if file !=tmpvar:
        system('mv '+tmpvar+' '+file)

def write_tnk(file,data):
    '''write TINKER xyz/arc
    '''
    file.write(" "+str(len(data))+"\n")
    for line in data:
       outline="%3s  %-3s%12.6f%12.6f%12.6f%6s%6s" % (line[0],line[1],float(line[2])\
                                     ,float(line[3]),float(line[4]),line[5],line[6])
       for i in range(7,len(line)):
          outline=outline+"%6s" % line[i]
       file.write(outline+"\n")
    return

def parsecmd():
    description="Build allotropic structures of Carbon, namely\
graphite hopg and armchair/zigzag carbon nanotubes.\n Output file can\
be saved in TINKER, XYZ, MOL2 or Gromacs GRO formats. Structures can also be periodic.\n"
    usage = "usage: %prog [options] output_file"
    #parse command line
    parser=OP(version='%prog 1.2',description=description, usage=usage)
    
    parser.add_option('-c','--credits',dest='credits',action='store_true',
                     default=False,help='display credits')
    parser.add_option('-s','--struct',dest='structure',default='none',
                      help='define structure: armcnt, zigzagcnt, hopg')
    parser.add_option('-p','--periodic',dest='pbc',action='store_true',
                     default=False, help='build periodic structure for Tinker')
    parser.add_option('-g','--geometry',dest='geometry',nargs=2,type='float',
                      help='define the geometry for the structure')
    parser.add_option('-f','--funct',dest='functionalization',default='none',
                      help='define functionalization: none, oh, cooh, coo- (only implemented for cnt)')
    parser.add_option('--xyz',dest='xyz',action='store_true',
                      help='write xyz file. This is FAST.')
    parser.add_option('--gro',dest='gro',action='store_true',
                      help='write gromacs gro file. This is FAST too ;).')
    parser.add_option('--mol2',dest='mol2',action='store_true',
                      help='write mol2 file.')
    (options, args) = parser.parse_args(argv[1:])
    
    #manage parse errors
    if options.credits: #display credits and quit
        credits="\n**********************************\n\
    Andrea Minoia, Martin Voegele\n\
    Contacts: minoiaa_at_gmail.com\
              http://chembytes.wikidot.com\
\n*********************************\n"
        print credits
        exit(0)

    if len(args)==0:   #arguments missing
        parser.exit(parser.print_help())
    
    if len(args)>1: #check if more than one argument (NOT OPTION) has been parsed
        parser.error('You have given me more than one argument '+str(args)+'... dunno what to do...\n')
    
    if lower(options.structure) != 'hopg' and lower(options.structure) !='armcnt'\
        and lower(options.structure) !='zigzagcnt':
        parser.error('Uknown structure: valid structures are hopg, armcnt and zigzagcnt')

    return options, args

def armcnt(n,l,ccbond,funct):
    ''' build armchair carbon nanotube
    '''
    atc=[]
    circ1=[]
    circ2=[]
    if funct=="oh":
        c1charge=-0.28
        c1acharge=0.24
        c2charge=0.01
    elif funct=="coo":
        c1charge=-0.34
        c1acharge=-0.09
        c2charge=0.03
    elif funct=="cooh":
        c1charge=-0.12
        c1acharge=-0.1
        c2charge=0.03
    else: 
        c1charge=-0.16
        c1acharge=-0.16
        c2charge=0.03
    cmcharge=0.0
    dx=ccbond*cos(120/2*(pi/180.0))
    dy=ccbond*sin(120/2*(pi/180.0))
    radius=(n*(2*dx+ccbond)+n*ccbond)/(2*pi)
    ycoord=+dy
    natoms=2*n
    #create circumferences
    for i in range(n):
        circ1.append(2*dx+ccbond)
        circ1.append(ccbond)
        circ2.append(ccbond)
        circ2.append(2*dx+ccbond)
    #adjust the circumferences
    circ1.insert(0,0.0)
    circ1.pop()
    circ2.insert(0,dx)
    circ2.pop()
    #Build CNT
    while ycoord>-l:
        ycoord-=dy
        arc=0.0
        # Assign suitable charges to the first and last two carbon rings
        if ycoord==0:
            ccharge_circ1a=c1acharge
            ccharge_circ1=c1charge
            ccharge_circ2a=c2charge
            ccharge_circ2=c2charge
        elif ycoord<-l+dy:
            ccharge_circ1a=c2charge
            ccharge_circ1=c2charge
            ccharge_circ2a=c1acharge
            ccharge_circ2=c1charge
        else:
            ccharge_circ1=cmcharge
            ccharge_circ2=cmcharge   
            ccharge_circ1a=cmcharge
            ccharge_circ2a=cmcharge 
        # Make coordinates
        for i in range(natoms):
            if modf(float(i)/2.0)[0]==0:
                newccharge=ccharge_circ1a
            else:
                newccharge=ccharge_circ1
            tmpcoords=['C']
            arc+=circ1[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
        ycoord-=dy
        arc=0.0
        for i in range(natoms):
            if modf(float(i)/2.0)[0]==0:
                newccharge=ccharge_circ2a
            else:
                newccharge=ccharge_circ2
            tmpcoords=['C']
            arc+=circ2[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
    
    pbc_l=abs(ycoord)+dy
    print '\n*******************************'
    print 'armchair CNT: n= ',n,' l (ang)= ',abs(ycoord)
    print 'periodicity (if apply) (ang)= ',pbc_l
    print 'diameter (ang): ',2*radius

    return atc,natoms,pbc_l,len(atc)

def zigzagcnt(n,l,ccbond,funct):
    ''' build zigzag carbon nanotube
    '''
    atc=[]
    circ1=[]
    circ2=[]

    if funct=="oh":
        c1charge=-0.10
        c1acharge=0.2
        c2charge=-0.04
    elif funct=="coo":
        c1charge=-0.43
        c1acharge=-0.33
        c2charge=0.09
    elif funct=="cooh":
        c1charge=-0.14
        c1acharge=-0.14
        c2charge=0.06
    else: 
        c1charge=-0.30
        c1acharge=-0.30
        c2charge=0.14
    cmcharge=0.0

    dy=ccbond*cos(120/2*(pi/180.0))
    dx=ccbond*sin(120/2*(pi/180.0))
    radius=(n*2*dx)/(2*pi)
    ycoord=+ccbond
    #create circumferences
    for i in range(n):
        circ1.append(2*dx)
    #adjust the circumferences
    
    circ1.pop()
    circ2=list(circ1)  #copy list!!! circ2=circ1 make duplicate and are both modified in the same way
    circ1.insert(0,0.0)
    circ2.insert(0,dx)

    #Build CNT
    while ycoord>-l:
        ycoord-=ccbond
        arc=0.0
	# Assign suitable charges to the first and last two carbon rings
        if ycoord==0:
            ccharge_circ1a=c1acharge
            ccharge_circ1=c1charge
            ccharge_circ2a=c2charge
            ccharge_circ2=c2charge
            ccharge_circ3a=cmcharge
            ccharge_circ3=cmcharge
            ccharge_circ4a=cmcharge
            ccharge_circ4=cmcharge
        elif ycoord<-l+dy:
            ccharge_circ1a=cmcharge
            ccharge_circ1=cmcharge
            ccharge_circ2a=cmcharge
            ccharge_circ2=cmcharge
            ccharge_circ3a=c2charge
            ccharge_circ3=c2charge
            ccharge_circ4a=c1acharge
            ccharge_circ4=c1charge
        else:
            ccharge_circ1=cmcharge
            ccharge_circ2=cmcharge  
            ccharge_circ3=cmcharge
            ccharge_circ4=cmcharge 
            ccharge_circ1a=cmcharge
            ccharge_circ2a=cmcharge 
            ccharge_circ3a=cmcharge
            ccharge_circ4a=cmcharge

        for i in range(n):
            if modf(float(i)/2.0)[0]==0:
                newccharge=ccharge_circ1a
            else:
                newccharge=ccharge_circ1
            tmpcoords=['C']
            arc+=circ1[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
        ycoord-=dy
        arc=0.0
        for i in range(n):
            if modf(float(i)/2.0)[0]==0:
                newccharge=ccharge_circ2a
            else:
                newccharge=ccharge_circ2
            tmpcoords=['C']
            arc+=circ2[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
        ycoord-=ccbond
        arc=0.0
        for i in range(n):
            if modf(float(i)/2.0)[0]==0:
                newccharge=ccharge_circ3a
            else:
                newccharge=ccharge_circ3
            tmpcoords=['C']
            arc+=circ2[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
        ycoord-=dy
        arc=0.0
        for i in range(n):
            if modf(float(i)/2.0)[0]==0:
                newccharge=ccharge_circ4a
            else:
                newccharge=ccharge_circ4
            tmpcoords=['C']
            arc+=circ1[i]
            theta=arc/radius
            tmpcoords.append(radius*cos(theta))
            tmpcoords.append(ycoord)
            tmpcoords.append(radius*sin(theta))
            tmpcoords.append('C.ar')
            tmpcoords.append(newccharge)
            atc.append(tmpcoords)
    pbc_l=abs(ycoord)+ccbond
    print '\n*******************************'
    print 'zigzag CNT: n= ',n,' l (ang)= ',abs(ycoord)
    print 'periodicity (if apply) (ang)= ',pbc_l
    print 'diameter (ang): ',2*radius

    return atc,n,pbc_l,len(atc)

def graphite(x,y,ccbond):
    ''' generate single square sheet of graphite HOPG
    '''
    atc=[]
    dx=ccbond*cos(120/2*(pi/180.0))
    dy=ccbond*sin(120/2*(pi/180.0))           
    ycoord=+dy
    xcoords1=[]
    xcoord=0.00
    xcoords1.append(xcoord)
    #build 1st row for X
    while xcoord<=x:
        xcoord+=ccbond+2*dx
        xcoords1.append(xcoord)
        xcoord+=ccbond
        xcoords1.append(xcoord)
    xcoords1.pop() #remove last element, i.e. the bond exceeding the size
    #build 2nd row for X
    xcoord=dx
    xcoords2=[]
    xcoords2.append(xcoord)
    while xcoord<=x+dx:
        xcoord+=ccbond
        xcoords2.append(xcoord)
        xcoord+=ccbond+2*dx
        xcoords2.append(xcoord)
    xcoords2.pop() #remove last element, i.e. the bond exceeding the size
    while ycoord>-y:
        ycoord-=dy
        for coord in xcoords1:
            tmpcoords=['C']
            tmpcoords.append(coord)
            tmpcoords.append(ycoord)
            tmpcoords.append(0.00)
            atc.append(tmpcoords)
        ycoord-=dy
        for coord in xcoords2:
            tmpcoords=['C']
            tmpcoords.append(coord)
            tmpcoords.append(ycoord)
            tmpcoords.append(0.00)
            atc.append(tmpcoords)
    print '\n*******************************'
    print 'HOPG graphite: a= ',atc[len(xcoords1)-1][1],' b= ',abs(ycoord)
    a_pbc=atc[len(xcoords1)-1][1]+ccbond
    b_pbc=abs(ycoord)+dy
    print 'Periodic (if apply) (ang): a= ',a_pbc, ' b= ',b_pbc

    return atc,len(xcoords1),a_pbc,b_pbc,len(atc)

def connect(coords,natx,pbcx,pbcy,nohcoords):
    '''build connectivity for graphite and nanotube:
    '''
    Ccov_r = 0.77  #covalent radius carbon
    Hcov_r = 0.32  #covalent radius Hydrogen
    Ocov_r = 0.66  #covalent radius Oxygen

    btollcc = (2*Ccov_r)*15/100  #bond tollerance of 15%
    btollch = (Hcov_r+Ccov_r)*5/100  #bond tollerance of 5% 
    btolloh = (Hcov_r+Ocov_r)*5/100  #bond tollerance of 5% 
    btollco = (Ccov_r+Ocov_r)*5/100  #bond tollerance of 5% 
    bondcc=[2*Ccov_r-btollcc,2*Ccov_r+btollcc]
    bondch=[(Hcov_r+Ccov_r)-btollch,(Hcov_r+Ccov_r)+btollch]
    bondco=[(Ocov_r+Ccov_r)-btollco,(Ocov_r+Ccov_r)+btollco]
    bondoh=[(Hcov_r+Ocov_r)-btolloh,(Hcov_r+Ocov_r)+btolloh]
    connect=zeros((len(coords),3),int) #init connectivity matrix
    bondlist=[]
    bondnumber=0
    #find connectivity, based on distance
    for i in xrange(len(coords)):
        for j in xrange(i+1,i+2*natx):
            if j<nohcoords:
                at1=[coords[i][1],coords[i][2],coords[i][3]]
                at2=[coords[j][1],coords[j][2],coords[j][3]]
                bond=getdist(at1,at2)
                if bond >= bondcc[0] and bond < bondcc[1]:
                    if i < j:
                        bondnumber+=1
                        tmpbond=[bondnumber,i+1,j+1,'ar']
                        bondlist.append(tmpbond)
                    for k in range(3):
                        if connect[i][k]== 0:
                            connect[i][k]=j+1   #index run from zero, not 1
                            break
                        else:
                            pass
                    for k in range(3):
                        if connect[j][k]== 0:
                            connect[j][k]=i+1   #index run from zero, not 1
                            break
                        else:
                            pass
        if len(coords)!=nohcoords: #there are hydrogens
            for j in xrange(nohcoords,len(coords)):
                at1=[coords[i][1],coords[i][2],coords[i][3]]
                at2=[coords[j][1],coords[j][2],coords[j][3]]
                bond=getdist(at1,at2)
                if ( ( bond >= bondch[0] and bond < bondch[1] ) or ( bond >= bondco[0] and bond < bondco[1] ) or ( bond >= bondoh[0] and bond < bondoh[1] ) or ( coords[i][0]=='C' and coords[j][0]=='C' and bond >= bondcc[0] and bond < bondcc[1] ) ) :
                    if i < j:
                        bondnumber+=1
                        tmpbond=[bondnumber,i+1,j+1,'1']
                        bondlist.append(tmpbond)
                    for k in range(3):
                        if connect[i][k]== 0:
                            connect[i][k]=j+1   #index run from zero, not 1
                            break
                        else:
                            pass
                    for k in range(3):
                        if connect[j][k]== 0:
                            connect[j][k]=i+1   #index run from zero, not 1
                            break
                        else:
                            pass
    #make periodic (TINKER) the structure                
    if pbcx:
        #cicle periodicity along X
        for i in xrange(len(coords)):
            for j in xrange(i,len(coords)):
                if coords[i][1]==0 and j-i==natx-1:
                    for k in range(3):
                        if connect[i][k]== 0:
                            connect[i][k]=j+1   #index run from zero, not 1
                            break
                        else:
                            pass
                    for k in range(3):
                        if connect[j][k]== 0:
                            connect[j][k]=i+1   #index run from zero, not 1
                            break
                        else:
                            pass
    if pbcy:
        i=0
        #cicle periodicity along Y
        while coords[i][2]==0.0:
            j=len(coords)-natx+i
            for k in range(3):
                if connect[i][k]== 0:
                    connect[i][k]=j+1   #index run from zero, not 1
                    break
                else:
                    pass
            for k in range(3):
                if connect[j][k]== 0:
                    connect[j][k]=i+1   #index run from zero, not 1
                    break
                else:
                    pass
            i+=1   
#    print connect, bondlist 
    return connect, bondlist
            
def data4tnk(coords,conn):
    ''' Prepare all the data to be written in
    tinker format
    '''
    data=[]
    for i in range(len(coords)):
        tmp=[i+1] #progressive index
        for j in range(4): #coordinates
            tmp.append(coords[i][j])
        tmp.append('0') #atom type
        for j in range(3): #connectivity
            if conn[i][j]==0:
                pass
            else:
               tmp.append(conn[i][j])
        data.append(tmp)
    return data

def add_COO(coords,natx,structure,is_protonated):
    '''Add COO- groups and hydrogens to nonperiodic structures
    '''
    Hcov_r = 0.32 
    Ocov_r = 0.66
    Ccov_r = 0.77

    chbond=Hcov_r+Ccov_r # 1.087 
    cobond=Ccov_r+Ocov_r # 1.362
    ohbond=Ocov_r+Hcov_r # 0.974
    ccbond=2*Ccov_r


    if structure =="hopg":  
        print "COO functionalization is not yet implemented for option hopg."
        exit(1)
                    
    elif structure == 'zigzagcnt':  

        if is_protonated:
            h1charge=0.19
	    c0charge=0.72
	    o1charge=-0.58
	    o2charge=-0.63
            hocharge=0.46
        else:
            h1charge=0.39
	    c0charge=0.91
	    o1charge=-0.86
	    o2charge=-0.86

	for i in range(len(coords)-natx,len(coords)):
            if modf(float(i)/2.0)[0]==0: #  COO-/COOH is added to every second C 
                Ox1=cobond*cos(i*2*pi/natx-pi/4)
	        Oz1=cobond*sin(i*2*pi/natx-pi/4)
                Ox2=cobond*cos(i*2*pi/natx+pi/4)
	        Oz2=cobond*sin(i*2*pi/natx+pi/4)
                #HOx=ohbond*cos(i*2*pi/natx+pi/4)
                #HOz=ohbond*sin(i*2*pi/natx+pi/4)
                tmpcoords=['C']
                tmpcoords.append(coords[i][1])
                tmpcoords.append(coords[i][2]-ccbond)
                tmpcoords.append(coords[i][3])
                tmpcoords.append('C.2')
		tmpcoords.append(c0charge)
                coords.append(tmpcoords)
                tmpcoords=['O']
		tmpcoords.append(coords[i][1]+Ox1) 
		tmpcoords.append(coords[i][2]-ccbond)
		tmpcoords.append(coords[i][3]+Oz1)
		tmpcoords.append('O.co2')
		tmpcoords.append(o1charge)
		coords.append(tmpcoords)
		tmpcoords=['O']
		tmpcoords.append(coords[i][1]+Ox2) 
		tmpcoords.append(coords[i][2]-ccbond)
		tmpcoords.append(coords[i][3]+Oz2)
		tmpcoords.append('O.co2')
		tmpcoords.append(o2charge)
		coords.append(tmpcoords)
                if is_protonated:
                    tmpcoords=['H']
		    tmpcoords.append(coords[i][1]+Ox2) 
		    tmpcoords.append(coords[i][2]-ccbond-ohbond)
		    tmpcoords.append(coords[i][3]+Oz2)
		    tmpcoords.append('H')
		    tmpcoords.append(hocharge)
		    coords.append(tmpcoords)    
            else:
                tmpcoords=['H']
                tmpcoords.append(coords[i][1])
                tmpcoords.append(coords[i][2]-chbond)
                tmpcoords.append(coords[i][3])
                tmpcoords.append('H')
		tmpcoords.append(h1charge)
                coords.append(tmpcoords)
        for i in range(natx):
            if modf(float(i)/2.0)[0]==0: #  COO- is added to every second C 
                Ox1=cobond*cos(i*2*pi/natx-pi/4)
	        Oz1=cobond*sin(i*2*pi/natx-pi/4)
                Ox2=cobond*cos(i*2*pi/natx+pi/4)
	        Oz2=cobond*sin(i*2*pi/natx+pi/4)
                #HOx=ohbond*cos(i*2*pi/natx+pi/4)
                #HOz=ohbond*sin(i*2*pi/natx+pi/4)
                tmpcoords=['C']
                tmpcoords.append(coords[i][1])
                tmpcoords.append(+ccbond)
                tmpcoords.append(coords[i][3])
                tmpcoords.append('C.2')
		tmpcoords.append(c0charge)
                coords.append(tmpcoords)
                tmpcoords=['O']
		tmpcoords.append(coords[i][1]+Ox1) 
		tmpcoords.append(+ccbond)
		tmpcoords.append(coords[i][3]+Oz1)
		tmpcoords.append('O.co2')
		tmpcoords.append(o1charge)
		coords.append(tmpcoords)
		tmpcoords=['O']
		tmpcoords.append(coords[i][1]+Ox2) 
		tmpcoords.append(+ccbond)
		tmpcoords.append(coords[i][3]+Oz2)
		tmpcoords.append('O.co2')
		tmpcoords.append(o2charge)
		coords.append(tmpcoords)
                if is_protonated:
                    tmpcoords=['H']
		    tmpcoords.append(coords[i][1]+Ox2) 
		    tmpcoords.append(coords[i][2]+ccbond+ohbond)
		    tmpcoords.append(coords[i][3]+Oz2)
		    tmpcoords.append('H')
		    tmpcoords.append(hocharge)
		    coords.append(tmpcoords)
            else:
                tmpcoords=['H']
                tmpcoords.append(coords[i][1])
                tmpcoords.append(+chbond)
                tmpcoords.append(coords[i][3])
                tmpcoords.append('H')
		tmpcoords.append(h1charge)
                coords.append(tmpcoords)


    elif structure == 'armcnt':

        if is_protonated:
            hocharge=0.44 # not assigned yet (all charges are diferent for coo and cooh)
            h1charge=0.17
	    c0charge=0.7
	    o1charge=-0.55
	    o2charge=-0.6
        else:
            h1charge=0.22
	    c0charge=0.83
	    o1charge=-0.84
	    o2charge=-0.84


	Hxz=chbond*cos(120/2*pi/180)
	Hy1=chbond*sin(120/2*pi/180)
	Cxz=ccbond*cos(120/2*pi/180)
	Cy1=ccbond*sin(120/2*pi/180)
	    
	for i in xrange(natx):  #upper border
	    Hx1=Hxz*cos(i*2*pi/natx)
	    Hz1=Hxz*sin(i*2*pi/natx)
	    Cx1=Cxz*cos(i*2*pi/natx)
	    Cz1=Cxz*sin(i*2*pi/natx)
	    Ox1=Cx1+cobond*cos(i*2*pi/natx - pi/4)
	    Oz1=Cz1+cobond*sin(i*2*pi/natx - pi/4)
	    Oy1=Cy1
	    Ox2=Cx1+cobond*cos(i*2*pi/natx + pi/4)
	    Oz2=Cz1+cobond*sin(i*2*pi/natx + pi/4)
	    Oy2=Cy1
            HOx=ohbond*cos(i*2*pi/natx+pi/4)
            HOz=ohbond*sin(i*2*pi/natx+pi/4)
            if modf(float(i)/2.0)[0]==0: #  COO- is added to every second C 
	#    if i==1:
		tmpcoords=['C']
		tmpcoords.append(coords[i][1]+Cx1) 
		tmpcoords.append(coords[i][2]+Cy1)
		tmpcoords.append(coords[i][3]+Cz1)
		tmpcoords.append('C.2')
		tmpcoords.append(c0charge)
		coords.append(tmpcoords)
		tmpcoords=['O']
		tmpcoords.append(coords[i][1]+Ox1) 
		tmpcoords.append(coords[i][2]+Oy1)
		tmpcoords.append(coords[i][3]+Oz1)
		tmpcoords.append('O.co2')
		tmpcoords.append(o1charge)
		coords.append(tmpcoords)
		tmpcoords=['O']
		tmpcoords.append(coords[i][1]+Ox2) 
		tmpcoords.append(coords[i][2]+Oy2)
		tmpcoords.append(coords[i][3]+Oz2)
		tmpcoords.append('O.co2')
		tmpcoords.append(o2charge)
		coords.append(tmpcoords)
                if is_protonated:
                    tmpcoords=['H']
		    tmpcoords.append(coords[i][1]+Ox2+HOx) 
		    tmpcoords.append(coords[i][2]+Oy2)
		    tmpcoords.append(coords[i][3]+Oz2+HOz)
		    tmpcoords.append('H')
		    tmpcoords.append(hocharge)
		    coords.append(tmpcoords)
	    else:
		tmpcoords=['H']
		tmpcoords.append(coords[i][1])        # +Hx1
		tmpcoords.append(coords[i][2]+chbond) # +Hy1
		tmpcoords.append(coords[i][3])        # +Hz1
		tmpcoords.append('H')
		tmpcoords.append(h1charge)
		coords.append(tmpcoords)

	if is_protonated:
            added_atoms=int(2.5*natx)
        else:
            added_atoms=int(2*natx) #  COO- is added to every second C 
#	    added_atoms=int(2+natx)	

        for i in xrange(len(coords)-added_atoms-natx,len(coords)-added_atoms): #bottom border
            Hx1=Hxz*cos(i*2*pi/natx)
            Hz1=Hxz*sin(i*2*pi/natx)
            Cx1=Cxz*cos(i*2*pi/natx)
            Cz1=Cxz*sin(i*2*pi/natx)
            Ox1=Cx1+cobond*cos(i*2*pi/natx - pi/4)
            Oz1=Cz1+cobond*sin(i*2*pi/natx - pi/4)
            Oy1=Cy1
            Ox2=Cx1+cobond*cos(i*2*pi/natx + pi/4)
            Oz2=Cz1+cobond*sin(i*2*pi/natx + pi/4)
            Oy2=Cy1
            HOx=ohbond*cos(i*2*pi/natx+pi/4)
            HOz=ohbond*sin(i*2*pi/natx+pi/4)
	    if modf(float(i)/2.0)[0]==0: #  COO- is added to every second C 
#	    if i==len(coords)-added_atoms-natx/2:
	        tmpcoords=['C']
                tmpcoords.append(coords[i][1]+Cx1) #update x even
	        tmpcoords.append(coords[i][2]-Cy1)
        	tmpcoords.append(coords[i][3]+Cz1)
	        tmpcoords.append('C.2')
	        tmpcoords.append(c0charge)
	        coords.append(tmpcoords)
	        tmpcoords=['O']
	        tmpcoords.append(coords[i][1]+Ox1) 
	        tmpcoords.append(coords[i][2]-Oy1)
	        tmpcoords.append(coords[i][3]+Oz1)
	        tmpcoords.append('O.co2')
	        tmpcoords.append(o1charge)
	        coords.append(tmpcoords)
	        tmpcoords=['O']
	        tmpcoords.append(coords[i][1]+Ox2) 
	        tmpcoords.append(coords[i][2]-Oy2)
	        tmpcoords.append(coords[i][3]+Oz2)
	        tmpcoords.append('O.co2')
	        tmpcoords.append(o2charge)
	        coords.append(tmpcoords)
                if is_protonated:
                    tmpcoords=['H']
		    tmpcoords.append(coords[i][1]+Ox2+HOx) 
		    tmpcoords.append(coords[i][2]-Oy2)
		    tmpcoords.append(coords[i][3]+Oz2+HOz)
		    tmpcoords.append('H')
		    tmpcoords.append(hocharge)
		    coords.append(tmpcoords)
	    else:
	        tmpcoords=['H']
	        tmpcoords.append(coords[i][1])        # +Hx1
	        tmpcoords.append(coords[i][2]-chbond) # -Hy1
	        tmpcoords.append(coords[i][3])        # +Hz1
	        tmpcoords.append('H')
	        tmpcoords.append(h1charge)
	        coords.append(tmpcoords)
           

def add_H(coords,natx,structure,funct_OH):
    '''
        Add hydrogens to nonperiodic structures
    '''

    Hcov_r = 0.32 
    Ocov_r = 0.66
    Ccov_r = 0.77
    
    chbond=Hcov_r+Ccov_r # 1.087 
    cobond=Ccov_r+Ocov_r # 1.362
    ohbond=Ocov_r+Hcov_r # 0.974

    Hxz=chbond*cos(120/2*pi/180)
    Hy1=chbond*sin(120/2*pi/180)
    Oxz=cobond*cos(120/2*pi/180)
    Oy1=cobond*sin(120/2*pi/180)

    if structure =="hopg" or structure=="armcnt": #saturate hopg in y direction or armchair cnt

        if funct_OH:
            h1charge=0.18
            o1charge=-0.53
            h2charge=0.37
        else:
            h1charge=0.13


        for i in xrange(natx):  #upper border
            Hx1=Hxz*cos(i*2*pi/natx)
            Hz1=Hxz*sin(i*2*pi/natx)
            Ox1=Oxz*cos(i*2*pi/natx)
            Oz1=Oxz*sin(i*2*pi/natx)
            Hx2=Ox1+ohbond*cos(i*2*pi/natx)
            Hz2=Oz1+ohbond*sin(i*2*pi/natx)
            Hy2=Oy1
            if structure=='armcnt': 

		if modf(float(i)/2.0)[0]==0 and funct_OH:
		    tmpcoords=['O']
		    tmpcoords.append(coords[i][1]+Ox1) 
		    tmpcoords.append(coords[i][2]+Oy1)
		    tmpcoords.append(coords[i][3]+Oz1)
		    tmpcoords.append('O.3')
		    tmpcoords.append(o1charge)
		    coords.append(tmpcoords)
		    tmpcoords=['H']
		    tmpcoords.append(coords[i][1]+Hx2) 
		    tmpcoords.append(coords[i][2]+Hy2)
		    tmpcoords.append(coords[i][3]+Hz2)
		    tmpcoords.append('H')
		    tmpcoords.append(h2charge)
		    coords.append(tmpcoords)
		else:
		    tmpcoords=['H']
		    tmpcoords.append(coords[i][1]+Hx1) 
		    tmpcoords.append(coords[i][2]+Hy1)
		    tmpcoords.append(coords[i][3]+Hz1)
		    tmpcoords.append('H')
		    tmpcoords.append(h1charge)
		    coords.append(tmpcoords)

            if structure=='hopg':   #set y and z for hopg
                tmpcoords=['H']
                tmpcoords.append(coords[i][1]+Hx1) #update x even
                tmpcoords.append(coords[i][2]+Hy1)
                tmpcoords.append(0.0)
	
	if funct_OH:
	    added_atoms=int(1.5*natx)
	else:
	    added_atoms=natx

        for i in xrange(len(coords)-added_atoms-natx,len(coords)-added_atoms): #bottom unsaturated border
            Hx1=Hxz*cos(i*2*pi/natx)
            Hz1=Hxz*sin(i*2*pi/natx)
            Ox1=Oxz*cos(i*2*pi/natx)
            Oz1=Oxz*sin(i*2*pi/natx)
            Hx2=Ox1+ohbond*cos(i*2*pi/natx)
            Hz2=Oz1+ohbond*sin(i*2*pi/natx)
            if structure=='armcnt':   #set y and z for hopg
		    if modf(float(i)/2.0)[0]==0 and funct_OH:
		        tmpcoords=['O']
		        tmpcoords.append(coords[i][1]+Ox1) #update x even
		        tmpcoords.append(coords[i][2]-Oy1)
		        tmpcoords.append(coords[i][3]+Oz1)
		        tmpcoords.append('O.3')
		        tmpcoords.append(o1charge)
		        coords.append(tmpcoords)
		        tmpcoords=['H']
		        tmpcoords.append(coords[i][1]+Hx2) #update x even
		        tmpcoords.append(coords[i][2]-Oy1)
		        tmpcoords.append(coords[i][3]+Hz2)
		        tmpcoords.append('H')
		        tmpcoords.append(h2charge)
		        coords.append(tmpcoords)
		    else:
		        tmpcoords=['H']
		        tmpcoords.append(coords[i][1]+Hx1) #update x odd
		        tmpcoords.append(coords[i][2]-Hy1)
		        tmpcoords.append(coords[i][3]+Hz1)
		        tmpcoords.append('H')
		        tmpcoords.append(h1charge)
		        coords.append(tmpcoords)

            if structure=='hopg':   #set y and z for hopg
                tmpcoords=['H']
                tmpcoords.append(coords[i][1]+Hx1) #update x even
                tmpcoords.append(coords[i][2]-Hy1)
                tmpcoords.append(0.0)
            
    if structure =="hopg":  #saturate hopg in x direction
        for i in xrange(len(coords)):
            if coords[i][0] == 'C':
                tmpcoords=['H']
                if coords[i][1]==0.0:
                    tmpcoords.append(-chbond)
                    tmpcoords.append(coords[i][2])
                    tmpcoords.append(0.0)
                    coords.append(tmpcoords)
                elif coords[i][1]==coords[natx-1][1]:
                    tmpcoords.append(coords[i][1]+chbond)
                    tmpcoords.append(coords[i][2])
                    tmpcoords.append(0.0)
                    coords.append(tmpcoords)
                    
    if structure == 'zigzagcnt':   #saturate zigzagcnt

        if funct_OH:
            h1charge=0.31
            o1charge=-0.40
            h2charge=0.07
        else:
            h1charge=0.16

        for i in range(len(coords)-natx,len(coords)):
            Hx=ohbond*cos(i*2*pi/natx)
            Hz=ohbond*sin(i*2*pi/natx)
            if modf(float(i)/2.0)[0]==0 and funct_OH:
		tmpcoords=['O']
		tmpcoords.append(coords[i][1]) #update x even
		tmpcoords.append(coords[i][2]-cobond)
		tmpcoords.append(coords[i][3])
		tmpcoords.append('O.3')
		tmpcoords.append(o1charge)
		coords.append(tmpcoords)
		tmpcoords=['H']
		tmpcoords.append(coords[i][1]+Hx) #update x even
		tmpcoords.append(coords[i][2]-cobond)
		tmpcoords.append(coords[i][3]+Hz)
		tmpcoords.append('H')
		tmpcoords.append(h2charge)
		coords.append(tmpcoords)
	    else:
                tmpcoords=['H']
                tmpcoords.append(coords[i][1])
                tmpcoords.append(coords[i][2]-chbond)
                tmpcoords.append(coords[i][3])
                tmpcoords.append('H')
		tmpcoords.append(h1charge)
                coords.append(tmpcoords)

        for i in range(natx):
            Hx=ohbond*cos(i*2*pi/natx)
            Hz=ohbond*sin(i*2*pi/natx)
            if modf(float(i)/2.0)[0]==0 and funct_OH:
		tmpcoords=['O']
		tmpcoords.append(coords[i][1]) #update x even
		tmpcoords.append(+cobond)
		tmpcoords.append(coords[i][3])
		tmpcoords.append('O.3')
		tmpcoords.append(o1charge)
		coords.append(tmpcoords)
		tmpcoords=['H']
		tmpcoords.append(coords[i][1]+Hx) #update x even
		tmpcoords.append(+cobond)
		tmpcoords.append(coords[i][3]+Hz)
		tmpcoords.append('H')
		tmpcoords.append(h2charge)
		coords.append(tmpcoords)
	    else:
                tmpcoords=['H']
                tmpcoords.append(coords[i][1])
                tmpcoords.append(+chbond)
                tmpcoords.append(coords[i][3])
                tmpcoords.append('H')
		tmpcoords.append(h1charge)
                coords.append(tmpcoords)



def write_xyz(file,data):
    '''
    Write a xyz file.
        Input Variables:
            file: output file (type: file)
            data: list of lists. Each list contains:
                    1) atom name
                2,3,4) X-, Y- and Z-coordinates
        Variables:
            line: store each list of data (type: list)
            outline: string containing a single line to be written in file (type: string)
    '''
    file.write(" "+str(len(data))+"\nGenerated by YASC buildCstruct v1.1\n")
    for line in data:
       outline="%-3s%12.6f%12.6f%12.6f" % (line[0],float(line[1]),float(line[2])\
                                     ,float(line[3]))
       file.write(outline+"\n")
    return

def write_gro(file,data,pbc1="",pbc2=""):
    '''
    Write a gromacs gro file.
        Input Variables:
            file: output file (type: file)
            data: list of lists. Each list contains:
                    1) atom name
                2,3,4) X-, Y- and Z-coordinates
            pbc1/pbc2: periodic lengths 
        Variables:
            line: store each list of data (type: list)
            outline: string containing a single line to be written in file (type: string)
    '''
    file.write("Generated by YASC buildCstruct v1.1\n "+str(len(data))+"\n")
    for index, line in enumerate(data):
       outline="%5i%-5s%5s%5i%8.3f%8.3f%8.3f" % (1,"CNT1",line[0],index,float(line[1])/10.0,float(line[2]/10.0)\
                                     ,float(line[3]/10.0))
       file.write(outline+"\n")
    if pbc1 == "":
        outline="  10   10   10\n"
    elif pbc2 == "":
        outline="  10   "+str(float(pbc1)/10.0)+"   10\n"
    else:
        outline="   "+str(float(pbc1)/10.0)+"   "+str(float(pbc2)/10.0)+"   1\n"
    file.write(outline+"\n")
    return


def write_mol2(file,data,bondlist):
    '''
    Write a mol2 file.
        Input Variables:
            file: output file (type: file)
            data: list of lists. Each list contains:
                    1) atom name
                2,3,4) X-, Y- and Z-coordinates
		    5) sybyl atom type
		    6) charge
            pbc1/pbc2: periodic lengths 
        Variables:
            line: store each list of data (type: list)
            outline: string containing a single line to be written in file (type: string)
    '''
    file.write("@<TRIPOS>MOLECULE\nCNT\n "+str(len(data))+" "+str(len(bondlist))+" 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
    for index, line in enumerate(data):
       outline="%7i %5s %8.3f %8.3f %8.3f %7s %7i %7s %8.3f" % (index+1,line[0],float(line[1]),float(line[2]),float(line[3]),line[4],1,"CNT1",float(line[5]))
       file.write(outline+"\n")

    file.write("@<TRIPOS>BOND\n")
    for line in bondlist:
       outline="%7i %7i %7i %7s" % (line[0],line[1],line[2],line[3])
       file.write(outline+"\n")
    return




def main():
    '''
    #===============================================================================
    #                               MAIN MAIN MAIN MAIN
    #===============================================================================
    '''
    ccbond=1.3874 # 1.42  #C-C bond length
    
    (options,args)=parsecmd()
    pbcx=False
    pbcy=False

    ofile=args[0] #get output file to save structure

    funct=lower(options.functionalization)

    if lower(options.structure) == "hopg":
        coords,natx,pbc_a,pbc_b,nohcoords=graphite(float(options.geometry[0]),float(options.geometry[1]),ccbond)
        if options.pbc:
            pbcx=True
            pbcy=True
        else:
            add_H(coords,natx,'hopg')
    elif lower(options.structure) == "armcnt":
        coords,natx,pbc_l,nohcoords=armcnt(int(options.geometry[0]),float(options.geometry[1]),ccbond,funct)
        if options.pbc:
            pbcx=False
            pbcy=True
        else:
            if funct == "coo":
		add_COO(coords,natx,'armcnt',False)
            elif funct == "cooh":
		add_COO(coords,natx,'armcnt',True)
            elif funct == "oh":
	        add_H(coords,natx,'armcnt',True)        
	    else:
	        add_H(coords,natx,'armcnt',False)
    else:
        coords,natx,pbc_l,nohcoords=zigzagcnt(int(options.geometry[0]),float(options.geometry[1]),ccbond,funct)
        if options.pbc:
            pbcx=False
            pbcy=True
        else:
            if funct == "coo":
		add_COO(coords,natx,'zigzagcnt',False)
            elif funct == "cooh":
		add_COO(coords,natx,'zigzagcnt',True)
	    elif funct == "oh":
	        add_H(coords,natx,'zigzagcnt',True)
	    else:
	        add_H(coords,natx,'zigzagcnt',False)
    print 'Atoms: ',len(coords)
    print 'saving structure...'
    if (not options.xyz and not options.gro) or options.mol2:
        conn,bondlist=connect(coords,natx,pbcx,pbcy,nohcoords) #get connectivity
        tnkdata=data4tnk(coords,conn) #write tinker output file
    backup_file(ofile)
    print '*******************************'
    OUT=open(ofile,'w')
    if options.xyz:
        write_xyz(OUT,coords)
    elif options.mol2:
        write_mol2(OUT,coords,bondlist)
    elif options.gro:
        if options.pbc:
            if lower(options.structure) == "hopg":
                write_gro(OUT,coords,pbc_a, pbc_b)
            else:
                write_gro(OUT,coords,pbc_l)
        else:
                write_gro(OUT,coords)
    else:
        write_tnk(OUT,tnkdata)
   
    exit(0)
    
if __name__ == "__main__":
    main()
