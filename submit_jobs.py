import os

# Get all the paths for the jobs
cwd = os.path.dirname(os.path.realpath(__file__))
abs_path="/".join(cwd.split('/')[:-1])
path_to_jobs="{}/xq1_jobs_test".format(abs_path)
run_set="cloud_order_xq1_test_fix_dt_2"
job_list=next(os.walk('{}'.format(path_to_jobs)))[1]
job_list.sort()

# for each job we want to go into that directory, then make and run rebound
for i in range(len(job_list)):
    run_dir="{}/{}/{}/{}".format(path_to_jobs,job_list[i],run_set,job_list[i])
    print "cd {}".format(run_dir)
    print "make"
    print "./rebound"
    print "\n"
