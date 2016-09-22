### This program sends the SQL request to the query and run the particle search job.
### Input file: halo catalog (argv[1]), job # (argv[2])
### Updated: Sept 3, 2015

import os, re, sys, commands, time
import numpy as np

### It takes care of the halos at the edge of the sim box, as the sim is periodic.
def Get_Piece(x_min, x_max, bit):
    if x_min > 0.0 and x_max < box_size: x_piece = "%s between %f and %f" % (bit, x_min, x_max)
    else:
        if x_min < 0.0:
            x_upper = x_min + box_size
            x_piece = "(%s > %f or %s < %f)" % (bit, x_upper, bit, x_max)
        else:
            x_lower = x_max - box_size
            x_piece = "(%s < %f or %s > %f)" % (bit, x_lower, bit, x_min)
    return x_piece

### login name and password for cosmosim
user = "xxx"
pw = "xxx"
web = "https://www.cosmosim.org/uws/query"
mock_link = "/Users/rickyccy/Documents/Research_assemblybias/mock/"
link = "/Users/rickyccy/Documents/Research_assemblybias/mock/halo/"

### You have to store the halo files in the directory "halo/"
#os.chdir(link + "halo/")
#file = np.loadtxt(sys.argv[1], unpack = True, dtype = {'names': ('id','x','y','z','R_vir'), 'formats':('S15','float64','float64','float64','float64')})

os.chdir(link)

"""
file = np.loadtxt(sys.argv[1], unpack = True, usecols = (0,1,2,3), dtype = {'names': ('id','x','y','z'), 'formats':('S15','float64','float64','float64')})

box_size = 250.0   ### box size of Bolshoi, 250.0 Mpc/h
halo_id = file[0]
x = file[1]
y = file[2]
z = file[3]
#dr = 1.2 * file[4]       ### Only concern about particles within R_vir in x-y plane.
#dz = 5. * file[4]   ### Go to 5 R_vir in z direction
dr = 1.0
dz = 1.0

del file

job_num = int(sys.argv[2])
halo_num = 1000

halo_len = len(x)

min_num = job_num * halo_num
max_num = (job_num + 1) * halo_num

if (max_num > halo_len):
    max_num = halo_len
"""

min_num = 0
max_num = 2

job_id = [0] * max_num

for i in range(min_num, max_num):
    job_id[i] = str(i).zfill(3)


#min_num = 834

### The text file to write down the job id's for downloading the particle files.
outfile = open(mock_link + "particles/job_id.txt", "w")

#outfile = open(mock_link + "particles/" + sys.argv[1][:3] + "/" + sys.argv[1][:-4] + "_" + sys.argv[2] + ".txt", "w")
#outfile = open(mock_link + "particles/" + sys.argv[1][-9:-4] + "_" + sys.argv[2] + ".txt", "w")

for i in range(min_num, max_num - 1):      ### replace 1 by len(x)
    ### Range of x, y and z of the concerned particles.
    """
    x_min = x[i] - dr   ## dr[i]
    x_max = x[i] + dr
    y_min = y[i] - dr
    y_max = y[i] + dr
    z_min = z[i] - dz
    z_max = z[i] + dz
    
    x_piece = Get_Piece(x_min, x_max, "p.x")
    y_piece = Get_Piece(y_min, y_max, "p.y")
    z_piece = Get_Piece(z_min, z_max, "p.z")
    """
    ### Job id name
    ### Firstly limit the range of x, y and z, then search for particles with z within n * R_vir from the center of the halo, where n = 5 in this case.
    #id_name = commands.getstatusoutput("uws --host %s --user %s --password %s job new \ query=\"SELECT p.x,p.y,p.z FROM Bolshoi.Particles416 AS p where %s and %s and %s and pdist(%6.1f,p.x,p.y,p.z,%f,%f,%f) < %f\" \ table=\"%s\" queue=\"long\" --run | sed -n -e 's/^.*Job ID: //p'" % (web, user, pw, x_piece, y_piece, z_piece, box_size, x[i], y[i], z[i], dz, halo_id[i]))

    id_name = commands.getstatusoutput("uws --host %s --user %s --password %s job new \ query=\"SELECT p.x,p.y,p.z FROM Bolshoi.Particles416 AS p WHERE particleId between 416%s0000000 and 416%s0000000-1\" \ table=\"%s\" queue=\"long\" --run | sed -n -e 's/^.*Job ID: //p'" % (web, user, pw, job_id[i], job_id[i+1], job_id[i]))

    outfile.write(str(id_name[1]) + "\t" + job_id[i] + "\n")
    
    #id_name = commands.getstatusoutput("uws --host %s --user %s --password %s job new \ query=\"SELECT rockstarId, pId, upId, Mvir, Rvir, Rs, x, y, z, vx,  vy, vz, Vmax, Macc, Mpeak, Vacc, Vpeak FROM MDPL2.Rockstar WHERE snapnum=125 AND Mvir between 1e10 AND 1e16\" \ table=\"%s\" queue=\"long\" --run " % (web, user, pw, i))
    
    print str(i) + " queries out of %d have been sent" % max_num

    time.sleep(7)

    ### Sleep for 6 mins for the previous queries to run. Avoid timeout.
    #if (i % 8 == 7):
    #    time.sleep(330)

outfile.close()


"""

#os.system("uws --host https://www.cosmosim.org/uws/query --user rickychue --password ilkitten0204 job run 368821763134053")

#os.system("uws --host https://www.cosmosim.org/uws/query --user rickychue --password ilkitten0204 job results 368821763134053 csv")

os.system("uws --host https://www.cosmosim.org/uws/query --user rickychue --password ilkitten0204 job delete list")
"""