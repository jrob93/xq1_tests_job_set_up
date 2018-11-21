import os
import subprocess
import numpy
import fileinput
import datetime

#create all the new run directories
base_dir="cloud_runs_fix_dt_2_base"
rebound_base_dir="../rebound_order_grav_collapse_condor_base"
dir_create_path="../xq1_jobs_test"
dir_name="cloud_order_xq1_test_fix_dt_2"
data_store="../../../../xq1_test_completed_jobs/"

# run parameters
num_list=[60,80,100]
X_rot=[0.5,0.5,0.5]
f_inf=[1,1,1]
R_eq=[100e3,250e3,750e3]
seed=["../../seed_pos/seed1_norm_10000.txt","../../seed_pos/seed1_norm_10000.txt","../../seed_pos/seed1_norm_10000.txt"]
Ntot=1e4
r_boulder = [-1.0,-1.0,-1.0]

collision_resolve= 'reb_collision_resolve_merge_out'
reb_func="../../reb_funcs/reb_func.h"
rebound_src="../.."

#variables to calculate dt and num_out
N_s=1 #number of files to sample binary orbit
# T_bin=5*24*60*60 #estimate of binary orbit period, s
#T_bin=100*24*60*60 #estimate of binary orbit period, s
T_bin=10*24*60*60 #estimate of binary orbit period, s
v_col=30 #expected collision velocity, m/s
t_max=1e0*365*24*60*60
t_max_str='1e0*365*24*60*60'

fname="run_deets_{}_list.txt".format(dir_name)
if os.path.isfile(fname):
    f=open(fname,"a")
else:
    f=open(fname,"w")
    f.write("num\tseed\tR_eq(m)\tX\tf\tr_b(m)\tdate\n")
print "list all jobs"
print "num,seed,R_eq,X,f"
for n in range(len(num_list)):
    print num_list[n],seed[n],R_eq[n],X_rot[n],f_inf[n]
    dt=f_inf[n]*R_eq[n]/(v_col*(Ntot**(1.0/3.0)))*(2.0/3.0)
    num_out=int(T_bin/dt/N_s)
    print dt
    print "Save a total of {} timesteps".format(t_max/dt/num_out)

    f.write("{:03d}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(num_list[n],seed[n],R_eq[n],X_rot[n],f_inf[n],r_boulder[n],str(datetime.date.today())))

f.close()

print "create jobs"

for n in range(len(num_list)):

    num=num_list[n]
    #create the new run directories from the base
    new_dir="{:03d}_{}".format(num,dir_name)
    print new_dir
    print num,seed[n],R_eq[n],X_rot[n],f_inf[n]

    reb_dir="{}/{}".format(dir_create_path,new_dir)
    run_set_dir="{}/{}".format(reb_dir,dir_name)
    run_dir="{}/{:03}_{}".format(run_set_dir,num,dir_name)

    # add safety check to not overwrite existing things?
    # add warning if we are overwriting?
    subprocess.call(["rm","-rf","{}/{}".format(dir_create_path,new_dir)])
    subprocess.call(["cp","-rf",rebound_base_dir,"{}/{}".format(dir_create_path,new_dir)])
    subprocess.call(["mkdir","{}".format(run_set_dir)])
    subprocess.call(["mkdir","{}/seed_pos".format(reb_dir)])
    subprocess.call(["mkdir","{}".format(run_dir)])

    subprocess.call("cp {}/* {}".format(base_dir,run_dir),shell=True)
    print "cp -rf reb_funcs {}".format(base_dir,run_dir)
    subprocess.call("cp -rf reb_funcs {}".format(reb_dir),shell=True)
    subprocess.call("cp seed_pos/{} {}/seed_pos".format(seed[n].split('/')[-1],reb_dir),shell=True)


    #set X,f,Req,seed,output interval...
    keyword1="double X"
    replacement1="double X={};\n".format(X_rot[n])

    keyword2="double _f"
    replacement2="double _f={};\n".format(f_inf[n])

    keyword3="rp->R_eq="
    replacement3="rp->R_eq={};\n".format(R_eq[n])

    keyword4="char ins[64]="
    replacement4="char ins[64]=\"{}\";\n".format(seed[n])

    #calculate num_out here from dt
    dt=f_inf[n]*R_eq[n]/(v_col*(Ntot**(1.0/3.0)))*(2.0/3.0)

    #CHECK NUM_OUT
    num_out=int(T_bin/dt/N_s)

    print "num_out: ",num_out
    print "Save a total of {} timesteps".format(t_max/dt/num_out)

    keyword5="double _num_out="
    replacement5="double _num_out={};\n".format(num_out)

    keyword6="double r_boulder="
    replacement6="double r_boulder={};\n".format(r_boulder[n])

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

    #fix Makefile
    keyword1="../../rebound_order"
    replacement1=rebound_src

    for line in fileinput.input("{}/Makefile".format(run_dir), inplace=True):
        if keyword1 in line:
            line = line.replace(keyword1, replacement1)
        print line,

    num+=1
    print "\n"

# Get all the paths for the jobs
cwd = os.path.dirname(os.path.realpath(__file__))
abs_path="/".join(cwd.split('/')[:-1])
path_to_jobs="{}/xq1_jobs_test".format(abs_path)
run_set="cloud_order_xq1_test_fix_dt_2"
job_list=next(os.walk('{}'.format(path_to_jobs)))[1]
job_list.sort()

# for each job we want to go into that directory, then make rebound, create submit.bash, and submit the job
for i in range(len(job_list)):
    run_dir="{}/{}/{}/{}".format(path_to_jobs,job_list[i],run_set,job_list[i])

    # compile the rebound binary
    print "make rebound"
    process=subprocess.Popen("make",shell=True,cwd="{}".format(run_dir)).wait()

    # Create the submit file from template
    print "create submit.bash"
    with open ("submit.bash", "r") as myfile:
        data=myfile.readlines()
    for j in range(len(data)):
        if '#SBATCH -J' in data[j]:
            data[j]='#SBATCH -J jamiedemo{}'.format(job_list[i].split('_')[0])

    # Add lines to copy back data and rm dir
    data.append('srun cp -rv {} {}\n'.format(run_dir,data_store))
    data.append('srun cd\n')
    data.append('srun rm -rfv {}/{}\n'.format(path_to_jobs,job_list[i]))

    with open ("{}/submit.bash".format(run_dir), "w") as myfile:
        myfile.writelines(data)

    # Submit job here
    process=subprocess.Popen("sbatch submit.bash",shell=True,cwd="{}".format(run_dir)).wait()
