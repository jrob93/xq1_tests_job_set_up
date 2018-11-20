import os
import subprocess

# Get all the paths for the jobs
cwd = os.path.dirname(os.path.realpath(__file__))
abs_path="/".join(cwd.split('/')[:-1])
path_to_jobs="{}/xq1_jobs_test".format(abs_path)
run_set="cloud_order_xq1_test_fix_dt_2"
job_list=next(os.walk('{}'.format(path_to_jobs)))[1]
job_list.sort()

# for each job create the submit.bash file from the default file
for i in range(len(job_list)):
    run_dir="{}/{}/{}/{}".format(path_to_jobs,job_list[i],run_set,job_list[i])

    with open ("submit.bash", "r") as myfile:
        data=myfile.readlines()

    for j in range(len(data)):
        # print data[j]
        if '#SBATCH -J' in data[j]:
            data[j]='#SBATCH -J jamiedemo{}'.format(job_list[i].split('_')[0])


    with open ("{}/submit.bash".format(run_dir), "w") as myfile:
        myfile.writelines(data)

    print data
    # process=subprocess.Popen("make",shell=True,cwd="{}".format(run_dir)).wait()
    # break
