import os
import subprocess
import numpy
import fileinput

#create all the new run directories
base_dir="cloud_runs_fix_dt_2_base"
rebound_base_dir="../rebound_order_grav_collapse_condor_base"
dir_create_path="../xq1_jobs_test"
dir_name="cloud_order_xq1_test_fix_dt_2"

#different inital conditions
# X_rot=[0.5,0.75,1.0,1.25]
# f_inf=[1,3,10,30,100]
# R_eq=[100e3,250e3,750e3]
# seed=["../../seed_pos/seed0_norm_100000.txt","../../seed_pos/seed1_norm_100000.txt","../../seed_pos/seed2_norm_100000.txt","../../seed_pos/seed3_norm_100000.txt"]
# r_boulder = 10.0
# collision_resolve= 'reb_collision_resolve_merge_out'
# reb_func="../../reb_funcs/reb_func.h"

X_rot=[0.75]
f_inf=[10]
R_eq=[250e3]
seed=["../../seed_pos/seed0_norm_10000.txt","../../seed_pos/seed1_norm_10000.txt","../../seed_pos/seed2_norm_10000.txt","../../seed_pos/seed3_norm_10000.txt"]
r_boulder = [-1.0]
collision_resolve= 'reb_collision_resolve_merge_out'
reb_func="../../reb_funcs/reb_func.h"
rebound_src="../.."

#variables to calculate dt and num_out
N_s=1 #number of files to sample binary orbit
# T_bin=5*24*60*60 #estimate of binary orbit period, s
#T_bin=100*24*60*60 #estimate of binary orbit period, s
T_bin=10*24*60*60 #estimate of binary orbit period, s
v_col=30 #expected collision velocity, m/s
Ntot=1e5
t_max=1e0*365*24*60*60
t_max_str='1e0*365*24*60*60'

# use this list to skip certain run numbers
skip_runs=[]

f=open("run_deets_{}.txt".format(dir_name),"w")
num=0
print "list all jobs"
print "num,seed,R_eq,X,f"
f.write("num\tseed\tR_eq(m)\tX\tf\tr_b(m)\n")
for g in range(len(seed)):
    for h in range(len(R_eq)):
        for i in range(len(X_rot)):
            for j in range(len(f_inf)):
                for k in range(len(r_boulder)):
                    print num,seed[g],R_eq[h],X_rot[i],f_inf[j]
                    dt=f_inf[j]*R_eq[h]/(v_col*(Ntot**(1.0/3.0)))*(2.0/3.0)
                    num_out=int(T_bin/dt/N_s)
                    print dt
                    print "Save a total of {} timesteps".format(t_max/dt/num_out)

                    f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(num,seed[g],R_eq[h],X_rot[i],f_inf[j],r_boulder[k]))
                    num+=1

f.close()

# exit()

# test=raw_input("Are you sure you want to rewrite all run files at {}? ".format(dir_create_path))
test='y'

if test == "y":

    num=0
    print "create jobs"
    for g in range(len(seed)):
        #if num>10:
        #    break
        for h in range(len(R_eq)):
            #if num>10:
            #    break
            for i in range(len(X_rot)):
                #if num>10:
                #    break
                for j in range(len(f_inf)):
                    #if num>10:
                    #    break
                    for k in range(len(r_boulder)):


                        #create the new run directories from the base
                        new_dir="{:03d}_{}".format(num,dir_name)
                        print new_dir
                        print num,seed[g],R_eq[h],X_rot[i],f_inf[j]

                        if num in skip_runs:
                            print 'skip this'
                            num+=1
                            print "\n"
                            continue

                        reb_dir="{}/{}".format(dir_create_path,new_dir)
                        run_set_dir="{}/{}".format(reb_dir,dir_name)
                        run_dir="{}/{:03}_{}".format(run_set_dir,num,dir_name)

                        subprocess.call(["rm","-rf","{}/{}".format(dir_create_path,new_dir)])

                        subprocess.call(["cp","-rf",rebound_base_dir,"{}/{}".format(dir_create_path,new_dir)])

                        subprocess.call(["mkdir","{}".format(run_set_dir)])
                        subprocess.call(["mkdir","{}/seed_pos".format(reb_dir)])

                        subprocess.call(["mkdir","{}".format(run_dir)])

                        subprocess.call("cp {}/* {}".format(base_dir,run_dir),shell=True)
                        print "cp -rf reb_funcs {}".format(base_dir,run_dir)
                        subprocess.call("cp -rf reb_funcs {}".format(reb_dir),shell=True)
                        subprocess.call("cp seed_pos/{} {}/seed_pos".format(seed[g].split('/')[-1],reb_dir),shell=True)


                        #set X,f,Req,seed,output interval...
                        keyword1="double X"
                        replacement1="double X={};\n".format(X_rot[i])

                        keyword2="double _f"
                        replacement2="double _f={};\n".format(f_inf[j])

                        keyword3="rp->R_eq="
                        replacement3="rp->R_eq={};\n".format(R_eq[h])

                        keyword4="char ins[64]="
                        replacement4="char ins[64]=\"{}\";\n".format(seed[g])

                        #calculate num_out here from dt
                        dt=f_inf[j]*R_eq[h]/(v_col*(Ntot**(1.0/3.0)))*(2.0/3.0)

                        #CHECK NUM_OUT
                        num_out=int(T_bin/dt/N_s)

                        print "num_out: ",num_out
                        print "Save a total of {} timesteps".format(t_max/dt/num_out)

                        keyword5="double _num_out="
                        replacement5="double _num_out={};\n".format(num_out)

                        keyword6="double r_boulder="
                        replacement6="double r_boulder={};\n".format(r_boulder[k])

                        keyword7="r->collision_resolve ="
                        replacement7="r->collision_resolve ={};\n".format(collision_resolve)

                        keyword8="#include \"../reb_funcs/reb_func.h\""
                        replacement8="#include \"{}\"\n".format(reb_func)

                        keyword9="static double t_max="
                        replacement9="static double t_max={};\n".format(t_max_str)

                        for line in fileinput.input("{}/problem.c".format(run_dir), inplace=True):
                            if keyword1 in line:
                                line = line.replace(line, replacement1)
                            if keyword2 in line:
                                line = line.replace(line, replacement2)
                            if keyword3 in line:
                                line = line.replace(line, replacement3)
                            if keyword4 in line:
                                line = line.replace(line, replacement4)
                            if keyword5 in line:
                                line = line.replace(line, replacement5)
                            if keyword6 in line:
                                line = line.replace(line, replacement6)
                            if keyword7 in line:
                                line = line.replace(line, replacement7)
                            if keyword8 in line:
                                line = line.replace(line, replacement8)
                            if keyword9 in line:
                                line = line.replace(line, replacement9)

                            print line,

                        #f = open("{}/{}/problem.c".format(dir_create_path,new_dir), 'r')
                        #print f.read()
                        #f.close()

                        #fix Makefile
                        keyword1="../../rebound_order"
                        replacement1=rebound_src

                        for line in fileinput.input("{}/Makefile".format(run_dir), inplace=True):
                            if keyword1 in line:
                                line = line.replace(keyword1, replacement1)
                            print line,

                        num+=1
                        print "\n"

        #                 break
        #             break
        #         break
        #     break
        # break
