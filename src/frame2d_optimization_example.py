from frame2d.frame2d import Frame2D, SteelBeam, SteelColumn, PointLoad,\
     FixedSupport, XHingedSupport, LineLoad, YHingedSupport, Hinge

import math
import random
import numpy as np

from scipy.optimize import fmin_slsqp, differential_evolution, fmin_cobyla
#from pso import pso


from timeit import default_timer as timer

from xlwt import Workbook,easyxf

import pylab
#import GA
#import datetime
#import pyvolution


# profiles from lightest to heaviest
PROFILES = [
                'IPE 80', 'IPE 100', 'IPE 120', 'HE 100 AA', 'IPE 140',
                'HE 120 AA', 'IPE 160','HE 100 A', 'HE 140 AA', 'IPE 180',
                'HE 120 A', 'HE 100 B', 'IPE 200', 'HE 160 AA','HE 140 A',
                'IPE 220', 'HE 120 B', 'HE 180 AA', 'HE 160 A', 'IPE 240',
                'HE 100 C','HE 140 B', 'HE 200 AA', 'HE 180 A', 'IPE 270',
                'HE 120 C', 'HE 220 AA', 'HE 100 M','IPE 300', 'HE 200 A',
                'HE 160 B', 'HE 240 AA', 'HE 140 C', 'IPE 330', 'HE 220 A',
                'HE 180 B', 'HE 120 M', 'HE 260 AA', 'IPE 360', 'HE 160 C',
                'HE 240 A', 'HE 280 AA','HE 200 B', 'HE 140 M', 'IPE 400',
                'HE 260 A', 'HE 300 AA', 'HE 180 C', 'HE 220 B','HE 320 AA',
                'HE 160 M', 'HE 280 A', 'IPE 450', 'HE 340 AA', 'HE 200 C',
                'HE 240 B','HE 360 AA', 'HE 300 A', 'HE 180 M', 'IPE 500',
                'HE 400 AA', 'HE 260 B', 'HE 220 C','HE 320 A', 'HE 450 AA',
                'HE 200 M', 'HE 280 B', 'HE 340 A', 'IPE 550', 'HE 500 AA',
                'HE 360 A', 'HE 300 B', 'HE 220 M', 'HE 240 C', 'HE 550 AA',
                'IPE 600', 'HE 400 A','HE 320 B', 'HE 600 AA', 'HE 260 C',
                'HE 340 B', 'HE 650 AA', 'HE 450 A', 'HE 360 B','HE 280 C',
                'HE 700 AA', 'HE 500 A', 'HE 400 B', 'HE 240 M', 'HE 550 A',
                'HE 450 B','HE 800 AA', 'HE 260 M', 'HE 300 C', 'HE 600 A',
                'HE 320 C', 'HE 500 B', 'HE 280 M','HE 650 A', 'HE 900 AA',
                'HE 550 B', 'HE 700 A', 'HE 600 B', 'HE 1000 AA', 'HE 800 A',
                'HE 650 B', 'HE 300 M', 'HE 700 B', 'HE 320 M', 'HE 340 M',
                'HE 360 M', 'HE 900 A','HE 400 M', 'HE 800 B', 'HE 450 M',
                'HE 500 M', 'HE 1000 A', 'HE 550 M', 'HE 600 M','HE 900 B',
                'HE 650 M', 'HE 700 M', 'HE 1000 B', 'HE 800 M', 'HE 900 M',
                'HE 1000 M','IPE 750']

IPE_PROFILES = ['IPE 80', 'IPE 100', 'IPE 120', 'IPE 140', 'IPE 160',
                'IPE 180', 'IPE 200', 'IPE 220', 'IPE 240', 'IPE 270',
                'IPE 300', 'IPE 330', 'IPE 360', 'IPE 400', 'IPE 450',
                'IPE 500', 'IPE 550', 'IPE 600', 'IPE 750']
#                       A mm^2            Iy mm^4
IPE_A_I_matrix = [[764.3401836602552, 801376.1398006069],
                  [1032.3219599741, 1710119.3948450384],
                  [1321.0219599741001, 3177531.491612866],
                  [1642.6019599741, 5412237.4783326965],
                  [2009.130995059227, 8692922.280149484],
                  [2394.7309950592266, 13169581.740920702],
                  [2848.41065788307, 19431661.97530533],
                  [3337.0506578830696, 27718364.6843977],
                  [3911.6216529422964, 38916214.50755847],
                  [4594.5016529422965, 57897773.31813715],
                  [5381.201652942297, 83561026.91741039], 
                  [6260.623980236907, 117668927.27986754],
                  [7272.9239802369075, 162656174.3911142],
                  [8446.3576397669, 231283456.02658707],
                  [9882.077639766901, 337429141.0374633],
                  [11552.1576397669, 481985026.72959846],
                  [13441.60263153228, 671164648.8112286],
                  [15598.44263153228, 920833980.6682992],
                  [17060.07972311255, 1507046150.557715]]

HEA_PROFILES = ['HE 100 A','HE 120 A','HE 140 A','HE 160 A','HE 180 A',
                'HE 200 A','HE 220 A','HE 240 A','HE 260 A','HE 280 A',
                'HE 300 A','HE 320 A','HE 340 A','HE 360 A','HE 400 A',
                'HE 450 A','HE 500 A','HE 550 A','HE 600 A','HE 650 A',
                'HE 700 A','HE 800 A','HE 900 A','HE 1000 A']

#                       A mm^2            Iy mm^4
HEA_A_I_matrix =[[2120, 3492000],
                 [2530, 6062000],
                 [3140, 10330000],
                 [3880, 16730000],
                 [4530, 25100000],
                 [5380, 36920000],
                 [6430, 54100000],
                 [7680, 77630000],
                 [8680, 104500000],
                 [9730, 136700000],
                 [11250, 182600000],
                 [12440, 229300000],
                 [13350, 276900000],
                 [14280, 330900000],
                 [15900, 450700000],
                 [17800, 637200000],
                 [19750, 869700000],
                 [21180, 1119000000],
                 [22650, 1412000000],
                 [24160, 1752000000],
                 [26050, 2153000000],
                 [28580, 3034000000],
                 [32050, 4221000000],
                 [34680,5538000000]]


HEB_PROFILES = ['HE 100 B','HE 120 B','HE 140 B','HE 160 B','HE 180 B',
                'HE 200 B','HE 220 B','HE 240 B','HE 260 B','HE 280 B',
                'HE 300 B','HE 320 B','HE 340 B','HE 360 B','HE 400 B',
                'HE 450 B','HE 500 B','HE 550 B','HE 600 B','HE 650 B',
                'HE 700 B','HE 800 B','HE 900 B','HE 1000 B']



HEB_A_I_matrix =[[2600, 4495000],
                 [3400,	8644000],
                 [4300,	15090000],
                 [5430,	24920000],
                 [6530,	38310000],
                 [7810,	56960000],
                 [9100,	80910000],
                 [10600, 112600000],
                 [11840, 149200000],
                 [13140, 192700000],
                 [14910, 251700000],
                 [16130, 308200000],
                 [17090, 366600000],
                 [18060, 431900000],
                 [19780, 576800000],
                 [21800, 798900000],
                 [23860, 1072000000],
                 [25410, 1367000000],
                 [27000, 1710000000],
                 [28630, 2106000000],
                 [30640, 2569000000],
                 [33420, 3591000000],
                 [37130, 4941000000],
                 [40000, 6447000000]]





def create_frame_coord():
    """
    Creates frame
    :return: Frame2D object
    """
    # Create your frame here
    # You can add parameters to the function call so your frame is created
    # according to those parameters

    # Coordinates of members, loads and supports
    
    
    coord1 = [[0.0,0], [0.0,3.5]]
    coord2 = [[0,3.5], [0,7]]
    coord3 = [[0,7], [0,10.5]]
    
    coord4 = [[6.0,0], [6.0,3.5]]
    coord5 = [[6,3.5], [6,7]]
    coord6 = [[6,7], [6,10.5]]
    
    coord7 = [[12.0,0], [12.0,3.5]]
    coord8 = [[12,3.5], [12,7]]
    coord9 = [[12,7], [12,10.5]]
    
    coord10 = [[18.0,0], [18.0,3.5]]
    coord11 = [[18,3.5], [18,7]]
    coord12 = [[18,7], [18,10.5]]
    
    
    coord13 = [[0.0,3.5], [6,3.5]]
    coord14 = [[0,7], [6,7]]
    coord15 = [[0,10.5], [6,10.5]]
    
    coord16 = [[6.0,3.5], [12,3.5]]
    coord17 = [[6,7], [12,7]]
    coord18 = [[6,10.5], [12,10.5]]
    
    coord19 = [[12.0,3.5], [18,3.5]]
    coord20 = [[12,7], [18,7]]
    coord21 = [[12,10.5], [18,10.5]]
    
    
    
    supp_coord1 = [0.0,0]
    supp_coord2 = [6,0]
    supp_coord3 = [12,0]
    supp_coord4 = [18,0]
    

    # Loads and members
    col1 = SteelColumn(coord1)
    col2 = SteelColumn(coord2)
    col3 = SteelColumn(coord3)
    
    col4 = SteelColumn(coord4)
    col5 = SteelColumn(coord5)
    col6 = SteelColumn(coord6)
    
    col7 = SteelColumn(coord7)
    col8 = SteelColumn(coord8)
    col9 = SteelColumn(coord9)
    
    col10 = SteelColumn(coord10)
    col11 = SteelColumn(coord11)
    col12 = SteelColumn(coord12)
    
       
    beam1 = SteelBeam(coord13)
    beam2 = SteelBeam(coord14)
    beam3 = SteelBeam(coord15)
    
    beam4 = SteelBeam(coord16)
    beam5 = SteelBeam(coord17)
    beam6 = SteelBeam(coord18)
    
    beam7 = SteelBeam(coord19)
    beam8 = SteelBeam(coord20)
    beam9 = SteelBeam(coord21)
    
    
    
    load1 = LineLoad(beam1, [-50.1,-50.1], 'y')
    load2 = LineLoad(beam2, [-50.1,-50.1], 'y')
    load3 = LineLoad(beam3, [-50.1,-50.1], 'y')
    load4 = LineLoad(beam4, [-50.1,-50.1], 'y')
    load5 = LineLoad(beam5, [-50.1,-50.1], 'y')
    load6 = LineLoad(beam6, [-50.1,-50.1], 'y')
    load7 = LineLoad(beam7, [-50.1,-50.1], 'y')
    load8 = LineLoad(beam8, [-50.1,-50.1], 'y')
    load9 = LineLoad(beam9, [-50.1,-50.1], 'y')
            
    load10 = PointLoad(coord1[1], [22.05, 0,0])
    load11 = PointLoad(coord2[1], [22.05, 0,0])
    load12 = PointLoad(coord3[1], [22.05, 0,0])
    
    load13 = PointLoad(coord10[1], [22.05, 0,0])
    load14 = PointLoad(coord11[1], [22.05, 0,0])
    load15 = PointLoad(coord12[1], [22.05, 0,0])
    
    

    # Create empty frame 'envelope'
    frame = Frame2D(num_elements=4)

    # Add members
    frame.add(col1)
    frame.add(col2)
    frame.add(col3)
    
    frame.add(col4)
    frame.add(col5)
    frame.add(col6)
    
    frame.add(col7)
    frame.add(col8)
    frame.add(col9)
    
    frame.add(col10)
    frame.add(col11)
    frame.add(col12)
    
    
    frame.add(beam1)
    frame.add(beam2)
    frame.add(beam3)
    
    frame.add(beam4)
    frame.add(beam5)
    frame.add(beam6)
    
    frame.add(beam7)
    frame.add(beam8)
    frame.add(beam9)
    

    # Add loads
    frame.add(load1)
    frame.add(load2)
    frame.add(load3)
    frame.add(load4)
    frame.add(load5)
    frame.add(load6)
    frame.add(load7)
    frame.add(load8)
    frame.add(load9)
    frame.add(load10)
    frame.add(load11)
    frame.add(load12)
    frame.add(load13)
    frame.add(load14)
    frame.add(load15)

    # Add supports
    frame.add(FixedSupport(supp_coord1, supp_id=1))
    frame.add(FixedSupport(supp_coord2, supp_id=1))
    frame.add(FixedSupport(supp_coord3, supp_id=1))
    frame.add(FixedSupport(supp_coord4, supp_id=1))
        

    #change column profiles to HEA
    
    for member in frame.members.values():
        if member.mtype=="column":
            member.profile = "HE 400 A"
            #member.profile = "IPE 100"
        else:
            member.profile = "HE 400 A"
            #member.profile = "IPE 330"

   
    # Generate frame
    frame.generate()
    # Return generated frame
   
    frame.plot()

    # Calculate result
    frame.calculate()
    #print(frame.nodal_forces[2])
    
    
    #frame.f.draw()
    frame.bmd(20)
    frame.plot_deflection()

    
    
    
    
    
    
    
    """
    coord1 = [[0.0,0], [0.0,1]]
    coord2 = [[0,1], [1,1]]
    coord3 = [[1,0], [1,1]]
    
    coord4 = [[0,1], [0,2]]
    coord5 = [[0,2], [1,2]]
    coord6 = [[1,1], [1,2]]
       
    supp_coord1 = [0.0,0]
    supp_coord2 = [1,0]

    # Loads and members
    col1 = SteelColumn(coord1)
    col2 = SteelColumn(coord3)
    beam1 = SteelBeam(coord2)
    
    
    col3 = SteelColumn(coord4)
    col4 = SteelColumn(coord6)
    beam2 = SteelBeam(coord5)
    
    
    load1 = PointLoad(coord1[1], [50, 0,0])
    load2 = LineLoad(beam1, [-1000,-1000], 'y')

    load3 = PointLoad(coord4[1], [50, 0,0])
    load4 = LineLoad(beam2, [-1500, -1500], 'y')
    

    # Create empty frame 'envelope'
    frame = Frame2D(num_elements=4)

    # Add members
    frame.add(col1)
    frame.add(col2)
    
    frame.add(col3)
    frame.add(col4)
    
    frame.add(beam1)
    frame.add(beam2)
    

    # Add loads
    frame.add(load1)
    frame.add(load2)
    frame.add(load3)
    frame.add(load4)

    # Add supports
    frame.add(FixedSupport(supp_coord1, supp_id=1))
    frame.add(FixedSupport(supp_coord2, supp_id=1))
        

    #change column profiles to HEA
    
    for member in frame.members.values():
        if member.mtype=="column":
            member.profile = "HE 400 A"
            #member.profile = "IPE 100"
        else:
            member.profile = "IPE 400"
            #member.profile = "IPE 330"

   
    # Generate frame
    frame.generate()
    # Return generated frame
   
    frame.plot()

    # Calculate result
    frame.calculate()
    #print(frame.nodal_forces[2])
    
    
    #frame.f.draw()
    frame.bmd(20)
    frame.plot_deflection()

    print ("weight =", frame.weight,"kg")
    """
       
   
   


    #for member in frame.members.values():
        #print ( "\n",member.profile, member.mtype,"h =",member.cross_section.h, "b=",member.cross_section.b, "tf=",
               #member.cross_section.tf,"tw=",member.cross_section.tw,"r=",member.cross_section.r)
    
        
    # Print nodal displacements and forces 
    
    #for i in range (len(frame.nodal_displacements)):
        #print (i,frame.nodal_displacements[i],frame.nodal_forces[i])
    
    
    sheet0 = book.add_sheet("Sheet0")
    bb = easyxf('border: bottom thick,top thin')
    bbl = easyxf('border: left thin, bottom thick,top thin')
    bbr = easyxf('border: right thin,bottom thick,top thin')
    bblr = easyxf('border: left thin, right thin, bottom thick,top thin')
    b = easyxf('border: bottom thin')
    bl = easyxf('border: left thin, bottom thin')
    br = easyxf('border: right thin, bottom thin')
    blr = easyxf('border: left thin, right thin, bottom thin')
    
    num=0
    sheet0.write(0, 0, "Member",style=bbl)
    sheet0.write(0, 1, "Element",style=bb) 
    sheet0.write(0, 2, "Disp X",style=bb) 
    sheet0.write(0, 3, "Disp Y",style=bb)
    sheet0.write(0, 4, "Disp Z",style=bb) 
    sheet0.write(0, 5, "N",style=bb) 
    sheet0.write(0, 6, "V",style=bb) 
    sheet0.write(0, 7, "M",style=bbr) 
    
    sheet0.write(0, 9, "Member",style=bbl)
    sheet0.write(0, 10, "Ratios",style=bbr)
    
        
    sheet0.write(0, 12, "Member",style=bbl)
    sheet0.write(0, 13, "Type",style=bb)
    sheet0.write(0, 14, "Profile",style=bb)
    sheet0.write(0, 15, "h",style=bb)
    sheet0.write(0, 16, "b",style=bb)
    sheet0.write(0, 17, "tf",style=bb)
    sheet0.write(0, 18, "tw",style=bb)
    sheet0.write(0, 19, "r",style=bbr)
    
    
    
    sheet0.write(0, 21, "Weight",style=bblr)
    #sheet0.write(3, 21, "Time",style=bblr)
    
    
     
    
    for i,member in enumerate (frame.members.values()):
        
        #print ("MEMBER:",i)
        #print ("Displacements member",i)
        
        row=sheet0.row(num+1)
        
        for key in (member.nodal_displacements):
            #print(key, "desp:",member.nodal_displacements[key], "forces:",member.nodal_forces[key])
            
            num+=1
            row=sheet0.row(num)
            row.write(0, i,style=bl)
            row.write(1, key,style=b)
            row.write(2, member.nodal_displacements[key][0],style=b)
            row.write(3, member.nodal_displacements[key][1],style=b)
            row.write(4, member.nodal_displacements[key][2],style=b)
            row.write(5, member.nodal_forces[key][0],style=b)
            row.write(6, member.nodal_forces[key][1],style=b)
            row.write(7, member.nodal_forces[key][2],style=br)
        
    
    num=0
    for i,member in enumerate (frame.members.values()):
        
        row=sheet0.row(num+1)
        for j in (member.r):
                       
            num+=1
            row=sheet0.row(num)
            row.write(9, i,style=bl)
            row.write(10, j,style=br)
            

        
    num=0        
    for i,member in enumerate (frame.members.values()):
        num+=1
        
        row=sheet0.row(num)
        row.write(12, i,style=bl)
        row.write(13, member.mtype,style=b)
        row.write(14, member.profile,style=b)
        row.write(15, member.cross_section.h,style=b)
        row.write(16, member.cross_section.b,style=b)
        row.write(17, member.cross_section.tf,style=b)
        row.write(18, member.cross_section.tw,style=b)
        row.write(19, member.cross_section.r,style=br)
     
    row=sheet0.row(1)
    row.write(21, frame.weight,style=blr)
   
    
    
    return frame
    


def create_frame_sym(bays,storeys,storey_height,bay_length):
    """
    Creates frame
    :return: Frame2D object
    """
    
    # Create empty frame 'envelope'
    frame = Frame2D(num_elements=4)    

    
    #Generate columns        

    # x-coordinate
    x = 0
    # y-coordinates 1: start, 2: end
    y1 = 0
    y2 = storey_height
    # number of columns
    num_cols = (bays+1)*storeys
    # create columns as FrameMember object
    for i in range(num_cols):
        frame.add(SteelColumn([[x, y1],[x, y2]]))  
                                                     
        # If the y2-coordinate equals frame's height, move the x-coordinate
        # by one bay length and initialize y-coordinates
        if y2 == storeys * storey_height:            
            x += bay_length
            y1 = 0
            y2 = storey_height
        else:
            y1 = y2
            y2 += storey_height
    
    
    #Generate beams
        
    # x-coordinate, 1:start, 2: end
    x1 = 0
    x2 = bay_length
    # y-coordinate
    y = storey_height
    # number of beams
    num_beams = bays*storeys
    # create beams as FrameMember object
    for i in range(num_beams):
        
        
        beam=SteelBeam([[x1,y],[x2,y]])
        frame.add(beam)
        
        """
        frame.add(LineLoad(beam, [-50.1,-50.1], 'y'))
        """
        
        # Generate lineloads
        
        if i==0:
            frame.add(LineLoad(beam, [-1000,-1000], 'y'))
            
        
        elif i==1:
            frame.add(LineLoad(beam, [-1500,-1500], 'y'))
            
       
        
        # If y-coordinate equals frame's height, initialize y-coordinate
        # and change x-coordinates
        if y == storeys * storey_height:
            y = storey_height
            x1 = x2
            x2 += bay_length  
        else:
            y += storey_height           
    
 
    

    # Add loads
    
       
    
    frame.add(PointLoad([0,1], [50, 0,0]))
    frame.add(PointLoad([0,2], [50, 0,0]))
   

    # Add supports
    x=0
    for i in range(bays+1):
        supp_coord = [x,0]
        x+=bay_length
        frame.add(FixedSupport(supp_coord, supp_id=1))
    

    #change column profiles to HEA
    
    for member in frame.members.values():
        if member.mtype=="column":
            member.profile = "HE 400 A"
            #member.profile = "IPE 400"
        else:
            member.profile = "IPE 400"
            #member.profile = "HE 400 B"

   
    # Generate frame
    frame.generate()
    # Return generated frame
   
    frame.plot()

    # Calculate result
    frame.calculate()
    #print(frame.nodal_forces[2])
    
    
    #frame.f.draw()
    frame.bmd(20)
    frame.plot_deflection()

    print ("weight =", frame.weight,"kg")
    
    
   
    sheet0 = book.add_sheet("Sheet0")
    bb = easyxf('border: bottom thick,top thin')
    bbl = easyxf('border: left thin, bottom thick,top thin')
    bbr = easyxf('border: right thin,bottom thick,top thin')
    bblr = easyxf('border: left thin, right thin, bottom thick,top thin')
    b = easyxf('border: bottom thin')
    bl = easyxf('border: left thin, bottom thin')
    br = easyxf('border: right thin, bottom thin')
    blr = easyxf('border: left thin, right thin, bottom thin')
    
    num=0
    sheet0.write(0, 0, "Member",style=bbl)
    sheet0.write(0, 1, "Element",style=bb) 
    sheet0.write(0, 2, "Disp X",style=bb) 
    sheet0.write(0, 3, "Disp Y",style=bb)
    sheet0.write(0, 4, "Disp Z",style=bb) 
    sheet0.write(0, 5, "N",style=bb) 
    sheet0.write(0, 6, "V",style=bb) 
    sheet0.write(0, 7, "M",style=bbr) 
    
    sheet0.write(0, 9, "Member",style=bbl)
    sheet0.write(0, 10, "Ratios",style=bbr)
    
        
    sheet0.write(0, 12, "Member",style=bbl)
    sheet0.write(0, 13, "Type",style=bb)
    sheet0.write(0, 14, "Profile",style=bb)
    sheet0.write(0, 15, "h",style=bb)
    sheet0.write(0, 16, "b",style=bb)
    sheet0.write(0, 17, "tf",style=bb)
    sheet0.write(0, 18, "tw",style=bb)
    sheet0.write(0, 19, "r",style=bbr)
    
    
    
    sheet0.write(0, 21, "Weight",style=bblr)
    #sheet0.write(3, 21, "Time",style=bblr)
    
    
     
    
    for i,member in enumerate (frame.members.values()):
        
        #print ("MEMBER:",i)
        #print ("Displacements member",i)
        
        row=sheet0.row(num+1)
        
        for key in (member.nodal_displacements):
            #print(key, "desp:",member.nodal_displacements[key], "forces:",member.nodal_forces[key])
            
            num+=1
            row=sheet0.row(num)
            row.write(0, i,style=bl)
            row.write(1, key,style=b)
            row.write(2, member.nodal_displacements[key][0],style=b)
            row.write(3, member.nodal_displacements[key][1],style=b)
            row.write(4, member.nodal_displacements[key][2],style=b)
            row.write(5, member.nodal_forces[key][0],style=b)
            row.write(6, member.nodal_forces[key][1],style=b)
            row.write(7, member.nodal_forces[key][2],style=br)
        
    
    num=0
    for i,member in enumerate (frame.members.values()):
        
        row=sheet0.row(num+1)
        for j in (member.r):
                       
            num+=1
            row=sheet0.row(num)
            row.write(9, i,style=bl)
            row.write(10, j,style=br)
            

        
    num=0        
    for i,member in enumerate (frame.members.values()):
        num+=1
        
        row=sheet0.row(num)
        row.write(12, i,style=bl)
        row.write(13, member.mtype,style=b)
        row.write(14, member.profile,style=b)
        row.write(15, member.cross_section.h,style=b)
        row.write(16, member.cross_section.b,style=b)
        row.write(17, member.cross_section.tf,style=b)
        row.write(18, member.cross_section.tw,style=b)
        row.write(19, member.cross_section.r,style=br)
     
    row=sheet0.row(1)
    row.write(21, frame.weight,style=blr)
    
    
    
    return frame
    

def cost_function(X):
    """
    Cost function used in optimization
    :param X: list of float values, e.g. member's dimensions
    :return cost_value: value of cost_function
    """
    # Implement cost function
    # Assign new optimized X values to the frame here
    # e.g. frame.members[1].h = X[1]
    
    
    #rho in kg/mm^3
   
    X = list(X)
    #print ("X :",X[0],X[5],X[10],X[15],X[20],X[25])
    for i, member in enumerate(frame.members.values()):
        #print (X[i])
        member.profile=X[i*5 : i*5 + 5]
        
        #print ("h dins cost:", i, member.cross_section.h,member.cross_section.tf)

        
        """
        member.h = X[i*5]
        member.cross_section.b = X[i*5+1]
        member.cross_section.tf = X[i*5+2]
        member.cross_section.tw = X[i*5+3]
        member.cross_section.r = X[i*5+4]
        """
        
        #print("abans de calcular:",i, member.length,member.profile, member.mtype, member.h,member.cross_section.h, member.cross_section.b,member.cross_section.tf)
      
                
    #frame.calculate()
   
    #for i, member in enumerate(frame.members.values()):
        #print ("MEMBRE ", i,member.cross_section.H,member.cross_section.h )   

  
  
    return frame.weight
    
    

   


def constraint_function(X):
    """
    Constraints used in the optimization
    :param X: list of float values, e.g. member's dimensions
    :return constraints: list of stress ratios
    """
    # Implement constraints for the optimization
    # Calculate new stress ratios with the new optimized dimensions X
    # e.g. stress_ratios = np.asarray(frame.r)
    # e.g. constraints = 1 - stress_ratios
    
    
        
    # Iterate through every member and change their height's
    
   
    X = list(X)
    
    for i, member in enumerate(frame.members.values()):
        #print (X[i])
        member.profile=X[i*5 : i*5 + 5]
      
    #print ("X constrain:",X)
    """
    for i, member in enumerate(frame.members.values()):
        #print ("X[i]",X[i])
        
        member.h = X[i*5]
        member.cross_section.b = X[i*5+1]
        member.cross_section.tf = X[i*5+2]
        member.cross_section.tw = X[i*5+3]
        #member.cross_section.r = X[i*5+4]
        member.cross_section.r = 1
    """
    frame.calculate()
    
    
    # Constraint values need to be >= 0
    # Stress indices must be less than 1
        
    #stress_indices = np.asarray(frame.r)
   
        #constf_values = 1 - stress_indices
    #constf_values = stress_indices - 1
             
    #print("constf_values", constf_values)
    
    #return constf_values.flatten()
    
  
    
    cons=[]
    
    for i,member in enumerate (frame.members.values()):
        height=0
        sway=0
        
        for j in range(len(member.nodal_coordinates)):
            if member.nodal_coordinates[j][1]>height:
                height=member.nodal_coordinates[j][1]
            
        for key in (member.nodal_displacements):
            if member.nodal_displacements[key][0]>sway:
                sway=frame.nodal_displacements[i][0]
                
        
        max_sway = height*1000 / 400 # SWAY LIMIT        
        cons1 = max_sway - sway
        cons.append(cons1)
        
        #print ("cons1", i, cons1)
        
     
 
    
      
    
    length=[]
     
    
    for i,member in enumerate (frame.members.values()):
        
        if member.mtype=="beam":
            
            max_y = 0
            length=(member.coordinates[1][0]-member.coordinates[0][0])*1000
            max_deflection =length / 300 #It can be different for each bay    
        
            for node in member.nodal_displacements.keys():
                y_val =-member.nodal_displacements[node][1]
                if y_val> max_y:
                    max_y = y_val
            cons2=max_deflection - max_y
            cons.append(cons2)
            
          
                                              
       
    # Cross-section, buckling and lateral-buckling strength
    for i,member in enumerate (frame.members.values()):
        
        for x in member.r:
            cons3=1-x
            cons.append(cons3)  
 
    
    
    for key in frame.members.keys():
        member = frame.members[key]
        b=member.cross_section.b
        h=member.cross_section.h
        tf=member.cross_section.tf
        tw=member.cross_section.tw
        
        #First upper boundary then lower boundary
        
        
        #print("b,h...",h,b,tf,tw)
        
        cons4=-b+1.0403*h+6.3 
        cons5=b-0.32*h-20 
        cons6=-tf+0.09*h+10.5
        cons7=tf-0.018*h-3.4
        cons8=-tw+0.0487*b+6.9 
        cons9=+tw-0.0174*b-1.7 
        
        #cons8=-tw+0.04*h+7.2 
        #cons9=+tw-0.014*h-2.4
        
        #print("const 4-9", cons4, cons5, cons6, cons7, cons8, cons9)
       
        cons.append(cons4)
        cons.append(cons5)
        cons.append(cons6)
        cons.append(cons7)
        cons.append(cons8)
        cons.append(cons9)
                
    return cons
    
    

def optimize_frame(c, debug=False):
    """
    Optimizes frame
    :param frame: Frame2D object
    :return:
    """
    # Implement function that optimizes the frame
    # Optimization functions usually call constraint_function and cost_function
    # in their own code
    
    
    #self.optimize_joints = joints
        # If no optimizer is given in the fuction call, use last used optimizer
        # Also checks if given optimizer is in the optimizer list and raises
        # error if not, otherwise changes frame's optimizer to that algorithm


    
    #"SLSQP":     
  
        
    #lb=[80,46,5.2,3.8,5]
    #ub=[1008,310,40,21,30]
    bounds= [(80,1008),(46,310),(5.2,40),(3.8,21),(5,30)]
    #bounds= [(80,1008),(46,310),(5,40),(4,21),(5,30)]
    maxiter=600
    
    # initial guess
    x0=[]
    bounds=bounds*len(frame.members)
    for i in bounds:
        x0.append(random.uniform(i[0],i[1]))
        
       
    if debug:
        print("x0: ",x0)
        iprint = 2
    else:
        iprint = 0
    
    start = timer()        
    out, fx, its, imode, smode = fmin_slsqp(cost_function,
                                                x0,
                                                f_ieqcons=constraint_function,
                                                bounds=bounds,
                                                iprint=iprint,
                                                iter=maxiter,
                                                full_output=True)
            
            
    constraint_function(out)
         
    end = timer()
    time= end - start   
        
    print("\n","Optimization",c,"\n","\n","Time elapsed: ", time, "s","\n","Number of iterations:",its,"\n","Continuous weight:",frame.weight,"\n")
    # print(constraint_function(out))
    
    r=0
            
    for i in constraint_function(out):
        if i<-0.05:
            r+=1
    if r>0:
        if c<10:
            c+=1
            print ("Coinstraints not satisfied", r)
            print (constraint_function(out))
            
            for i, member in enumerate(frame.members.values()):
                print(i, member.profile, member.h, member.cross_section.b, member.cross_section.tf, member.cross_section.tw,member.cross_section.r,"\n")
        
            optimize_frame(c)
        else:
            print ("Coinstraints not satisfied", r)
            print ("Not feasible solution")
            
            
      
    
    else:                                         
        for i, member in enumerate(frame.members.values()):
            print(i, member.profile, member.h, member.cross_section.h, member.cross_section.b, member.cross_section.tf, member.cross_section.tw,member.cross_section.r,"\n")
            
            
        print("\n","Optimized weight:",frame.weight)    
        
         
        #for i, member in enumerate(frame.members.values()):
            #print(i, member.profile, member.mtype, member.h, member.cross_section.b,  member.cross_section.H, member.cross_section.B)
        
        
       
        frame.calculate()
        
       
        
        frame.plot()
        #pylab.savefig('frame_continuous.png')
        
        frame.bmd(20)
        frame.plot_deflection(10)
       
        
    
        #SAVE RESULTS
        
        #book = Workbook()
                
        sheet1 = book.add_sheet("Sheet1")
        bb = easyxf('border: bottom thick,top thin')
        bbl = easyxf('border: left thin, bottom thick,top thin')
        bbr = easyxf('border: right thin,bottom thick,top thin')
        bblr = easyxf('border: left thin, right thin, bottom thick,top thin')
        b = easyxf('border: bottom thin')
        bl = easyxf('border: left thin, bottom thin')
        br = easyxf('border: right thin, bottom thin')
        blr = easyxf('border: left thin, right thin, bottom thin')
        
        num=0
        sheet1.write(0, 0, "Member",style=bbl)
        sheet1.write(0, 1, "Element",style=bb) 
        sheet1.write(0, 2, "Disp X",style=bb) 
        sheet1.write(0, 3, "Disp Y",style=bb)
        sheet1.write(0, 4, "Disp Z",style=bb) 
        sheet1.write(0, 5, "N",style=bb) 
        sheet1.write(0, 6, "V",style=bb) 
        sheet1.write(0, 7, "M",style=bbr) 
        
        sheet1.write(0, 9, "Member",style=bbl)
        sheet1.write(0, 10, "Ratios",style=bbr)
        
            
        sheet1.write(0, 12, "Member",style=bbl)
        sheet1.write(0, 13, "Type",style=bb)
        sheet1.write(0, 14, "Profile",style=bb)
        sheet1.write(0, 15, "h",style=bb)
        sheet1.write(0, 16, "b",style=bb)
        sheet1.write(0, 17, "tf",style=bb)
        sheet1.write(0, 18, "tw",style=bb)
        sheet1.write(0, 19, "r",style=bbr)
        
        
        
        sheet1.write(0, 21, "Weight",style=bblr)
        sheet1.write(3, 21, "Time",style=bblr)
        sheet1.write(6, 21, "Its",style=bblr)
        
        
         
        
        for i,member in enumerate (frame.members.values()):
            
            #print ("MEMBER:",i)
            #print ("Displacements member",i)
            
            row=sheet1.row(num+1)
            
            for key in (member.nodal_displacements):
                #print(key, "desp:",member.nodal_displacements[key], "forces:",member.nodal_forces[key])
                
                num+=1
                row=sheet1.row(num)
                row.write(0, i,style=bl)
                row.write(1, key,style=b)
                row.write(2, member.nodal_displacements[key][0],style=b)
                row.write(3, member.nodal_displacements[key][1],style=b)
                row.write(4, member.nodal_displacements[key][2],style=b)
                row.write(5, member.nodal_forces[key][0],style=b)
                row.write(6, member.nodal_forces[key][1],style=b)
                row.write(7, member.nodal_forces[key][2],style=br)
            
        
        num=0
        for i,member in enumerate (frame.members.values()):
            
            row=sheet1.row(num+1)
            for j in (member.r):
                           
                num+=1
                row=sheet1.row(num)
                row.write(9, i,style=bl)
                row.write(10, j,style=br)
                
    
            
        num=0        
        for i,member in enumerate (frame.members.values()):
            num+=1
            
            row=sheet1.row(num)
            row.write(12, i,style=bl)
            row.write(13, member.mtype,style=b)
            row.write(14, member.profile,style=b)
            row.write(15, member.cross_section.h,style=b)
            row.write(16, member.cross_section.b,style=b)
            row.write(17, member.cross_section.tf,style=b)
            row.write(18, member.cross_section.tw,style=b)
            row.write(19, member.cross_section.r,style=br)
         
        row=sheet1.row(1)
        row.write(21, frame.weight,style=blr)
        
        row=sheet1.row(4)
        row.write(21, time,style=blr)
        
        row=sheet1.row(7)
        row.write(21, its,style=blr)
        
       
        
        #sheet1.insert_bitmap('frame_continuous.jpg',0,20)
        
        
           
        #book.save("test.xls") 
       


def check_constrains_discrete():
    
    frame.calculate()    
    
    change_profile=[] #list of the position of the members that don't satisfy the constrains
    cons=[]
    c=0
        
    for i,member in enumerate (frame.members.values()):
        height=0
        sway=0
        
        for j in range(len(member.nodal_coordinates)):
            if member.nodal_coordinates[j][1]>height:
                height=member.nodal_coordinates[j][1]
            
        for key in (member.nodal_displacements):
            if member.nodal_displacements[key][0]>sway:
                sway=frame.nodal_displacements[i][0]
                
        
        max_sway = height*1000 / 400 # SWAY LIMIT        
        cons1 = max_sway - sway
        cons.append(cons1)
        
        #print ("cons1", i, cons1)
        
        if cons1 < -0.05:  
            change_profile.append(i)
            c+=1
            
   
    
  
    # Beam deflection
                     
    length=[]
     
    
    for i,member in enumerate (frame.members.values()):
        
        if member.mtype=="beam":
            
            max_y = 0
            length=(member.coordinates[1][0]-member.coordinates[0][0])*1000
            max_deflection =length / 300 #It can be different for each bay    
        
            for node in member.nodal_displacements.keys():
                y_val =-member.nodal_displacements[node][1]
                if y_val> max_y:
                    max_y = y_val
            cons2=max_deflection - max_y
            cons.append(cons2)
            
            if cons2 < -0.05:
                if i not in change_profile:
                    change_profile.append(i) 
                    c+=1
   
          
                                              
       
    # Cross-section, buckling and lateral-buckling strength
    for i,member in enumerate (frame.members.values()):
        
        for x in member.r:
            cons3=1-x
            cons.append(cons3)  
            #print ("const3",cons3)
            
            if cons3 < -0.05:
                if i not in change_profile:
                    change_profile.append(i) 
                    c+=1
                    
    #print("r",key,frame.r)
    #print("cons r",cons)

   
            
    return c,change_profile

    
            
     
 

def continuous_to_discrete():
    """
    Changes frame's member's continuous dimensions to discrete profiles
    :param frame: optimized Frame2D object
    """
    
    # Implement function that changes continuously optimized members to discrete profiles
    
    start = timer()
    
    h_opt=[]
    b_opt=[]
    tf_opt=[]
    tw_opt=[]
    
    
    r_opt=[]
    
    a_opt=[]
    i_opt=[]
    
    
    for member in frame.members.values():
        h_opt.append(member.h)
        b_opt.append(member.cross_section.b)
        tf_opt.append(member.cross_section.tf)
        tw_opt.append(member.cross_section.tw)
        
        
        r_opt.append(member.cross_section.r)
        
        #h=member.h
        #b=member.cross_section.b
        #tf=member.cross_section.tf
        #tw=member.cross_section.tw
        #r=member.cross_section.r
        
        #a=2.0*tf*b + (h - 2.0*tf)*tw + (4.0 - math.pi)*r**2.0
        #a_opt.append(a)
        #iner=1.0/12.0*(b*h**3.0 - (b - tw)*(h - 2.0*tf)**3.0) + 0.03*r**4 + 0.2146*r**2.0*(h - 2.0*tf - 0.4468*r)**2.0
        #i_opt.append(iner)
    
    dic= {} 
    
    for i,member in enumerate(frame.members.values()):
        
        distances=[]
        
        for j in catalogue_profile.keys():
            hj=catalogue_profile[j]["h"]
            bj=catalogue_profile[j]["b"]
            tfj=catalogue_profile[j]["t_f"]
            twj=catalogue_profile[j]["t_w"]
            rj=catalogue_profile[j]["r"]
            
            #aj=2.0*tfj*bj + (hj - 2.0*tfj)*twj + (4.0 - math.pi)*rj**2.0
            #dj=aj
            
            #ij=1.0/12.0*(bj*hj**3.0 - (bj - twj)*(hj - 2.0*tfj)**3.0) + 0.03*rj**4 + 0.2146*rj**2.0*(hj - 2.0*tfj - 0.4468*rj)**2.0
            
            #dj=np.sqrt(((aj-a_opt[i])/a_opt[i])**2+((ij-i_opt[i])/i_opt[i])**2)
            
                        
            
            
            
            dj=np.sqrt(((hj-h_opt[i])/h_opt[i])**2+((bj-b_opt[i])/b_opt[i])**2+((tfj-tf_opt[i])/tf_opt[i])**2+((twj-tw_opt[i])/tw_opt[i])**2)
            
            #dj=np.sqrt(((hj-h_opt[i])/h_opt[i])**2+((bj-b_opt[i])/b_opt[i])**2)
            
            
            distances.append((dj,j))
        
        distances.sort()
        
         
        
        dic[i]=distances
        
    print ("dict:",dic)
        
             
        
    for i,member in enumerate(frame.members.values()):
        member.profile=dic[i][0][1]
        
        
        #print(i, member.profile, member.h, member.cross_section.b, member.cross_section.tf, member.cross_section.tw)
    #print("Discrete weight:",frame.weight)
            
   
    frame.calculate()
    
    c,change_profile=check_constrains_discrete()
    
    end=timer()
    time= end - start 
    
    
    p=[]
    
    for member in frame.members.values():
        p.append(0)
    
    
    while c>0:
        
        print ("Coinstraints not satisfied")
        print ("Change profile member:",c,change_profile, p)
        change_profile1=change_profile
        
        
      
        
        for i, member in enumerate(frame.members.values()):
            if i in change_profile and i in change_profile1:
                p[i]+=1
                print(" i in change profile",i)
                member.profile=dic[i][int(p[i])][1] 
                
                print(i, member.profile)
                
                #print(member.h, member.cross_section.h,catalogue_profile[member.profile]["h"], member.cross_section.b, member.cross_section.tf, member.cross_section.tw)
                
                """
                member.cross_section.h= catalogue_profile[member.profile]["h"]
                member.cross_section.b= catalogue_profile[member.profile]["b"]
                member.cross_section.tf= catalogue_profile[member.profile]["t_f"]
                member.cross_section.tw= catalogue_profile[member.profile]["t_w"]
                member.cross_section.r= catalogue_profile[member.profile]["r"]
                member.cross_section.catalogue=True
                """
                
                c1,change_profile1=check_constrains_discrete()
                #print(frame.r)
                #print ("PES:",frame.weight, c1,change_profile1)
                if c1==0:
                    c=c1
                    print ("ja compleix")
                    break # vui que es pari el bucle
        
        
        c=c1
        change_profile=change_profile1
        
    
    
    
    
    
    
    
  
    frame.calculate()            
    end=timer()
    time= end - start   
       
    for i, member in enumerate(frame.members.values()):
        print(i, member.profile,member.cross_section.h, member.cross_section.b, member.cross_section.tf, member.cross_section.tw)
            
       
    print("\n","Time elapsed: ", time, "s","\n","Discrete weight:",frame.weight)
            
        
    
    frame.plot()
    frame.bmd(20)
    frame.plot_deflection(10)
    
        
        
        
    #SAVE RESULTS
    
    #book = Workbook()
    
    sheet2 = book.add_sheet("Sheet2")
    bb = easyxf('border: bottom thick,top thin')
    bbl = easyxf('border: left thin, bottom thick,top thin')
    bbr = easyxf('border: right thin,bottom thick,top thin')
    bblr = easyxf('border: left thin, right thin, bottom thick,top thin')
    b = easyxf('border: bottom thin')
    bl = easyxf('border: left thin, bottom thin')
    br = easyxf('border: right thin, bottom thin')
    blr = easyxf('border: left thin, right thin, bottom thin')
    
    num=0
    sheet2.write(0, 0, "Member",style=bbl)
    sheet2.write(0, 1, "Element",style=bb) 
    sheet2.write(0, 2, "Disp X",style=bb) 
    sheet2.write(0, 3, "Disp Y",style=bb)
    sheet2.write(0, 4, "Disp Z",style=bb) 
    sheet2.write(0, 5, "N",style=bb) 
    sheet2.write(0, 6, "V",style=bb) 
    sheet2.write(0, 7, "M",style=bbr) 
    
    sheet2.write(0, 9, "Member",style=bbl)
    sheet2.write(0, 10, "Ratios",style=bbr)
    
        
    sheet2.write(0, 12, "Member",style=bbl)
    sheet2.write(0, 13, "Type",style=bb)
    sheet2.write(0, 14, "Profile",style=bb)
    sheet2.write(0, 15, "h",style=bb)
    sheet2.write(0, 16, "b",style=bb)
    sheet2.write(0, 17, "tf",style=bb)
    sheet2.write(0, 18, "tw",style=bb)
    sheet2.write(0, 19, "r",style=bbr)
    
    
    
    sheet2.write(0, 21, "Weight",style=bblr)
    sheet2.write(3, 21, "Time",style=bblr)
    
    
     
    
    for i,member in enumerate (frame.members.values()):
        
        #print ("MEMBER:",i)
        #print ("Displacements member",i)
        
        row=sheet2.row(num+1)
        
        for key in (member.nodal_displacements):
            #print(key, "desp:",member.nodal_displacements[key], "forces:",member.nodal_forces[key])
            
            num+=1
            row=sheet2.row(num)
            row.write(0, i,style=bl)
            row.write(1, key,style=b)
            row.write(2, member.nodal_displacements[key][0],style=b)
            row.write(3, member.nodal_displacements[key][1],style=b)
            row.write(4, member.nodal_displacements[key][2],style=b)
            row.write(5, member.nodal_forces[key][0],style=b)
            row.write(6, member.nodal_forces[key][1],style=b)
            row.write(7, member.nodal_forces[key][2],style=br)
        
    
    num=0
    for i,member in enumerate (frame.members.values()):
        
        row=sheet2.row(num+1)
        for j in (member.r):
                       
            num+=1
            row=sheet2.row(num)
            row.write(9, i,style=bl)
            row.write(10, j,style=br)
            

        
    num=0        
    for i,member in enumerate (frame.members.values()):
        num+=1
        
        row=sheet2.row(num)
        row.write(12, i,style=bl)
        row.write(13, member.mtype,style=b)
        row.write(14, member.profile,style=b)
        row.write(15, member.cross_section.h,style=b)
        row.write(16, member.cross_section.b,style=b)
        row.write(17, member.cross_section.tf,style=b)
        row.write(18, member.cross_section.tw,style=b)
        row.write(19, member.cross_section.r,style=br)
     
    row=sheet2.row(1)
    row.write(21, frame.weight,style=blr)
    row=sheet2.row(4)
    row.write(21, time,style=blr)
   
    #book.save("test.xls")
     
        
     

if __name__ == '__main__':
    # Implement here function calls to optimize the frame
    """
    e.g.
    frame = create_frame()
    optimize_frame(frame)
    continuous_to_discrete(frame)
    frame.plot()
    """
    book = Workbook()
    
    # create_frame_coord() entering all the coordinates
   
    frame=create_frame_coord()
    
    
    # create_frame_sym simetric frame enterin number of bays, height
   
    
    
    bays=1
    storeys=2
    storey_height=1
    bay_length=1
    #frame=create_frame_sym(bays,storeys,storey_height,bay_length) 
    #frame=create_frame_coord() 
    book.save("test.xls")
   
    
  
    print ("\n","CONTINUOUS PROBLEM")
    
    c=1
    optimize_frame(c, debug=True)
    book.save("test.xls")
    
    
    
    print ("\n","DISCRETE PROBLEM","\n")
    
    continuous_to_discrete()
    
    book.save("test.xls")
    
     
        
    #for member in frame.members.values():
         #print (member.profile,"h =",member.cross_section.h, "b=",member.cross_section.b, "tf=",
               #member.cross_section.tf,"tw=",member.cross_section.tw,"r=",member.cross_section.r)
    
    
     # Print nodal displacements and forces 
    """
    for i in range (len(frame.nodal_displacements)):
        print (i,frame.nodal_displacements[i],frame.nodal_forces[i])
    """
    
    
        

    
   
    