import subprocess

# process1=subprocess.Popen("python create_runs.py",shell=True).wait()
process1=subprocess.Popen("python create_a_run.py",shell=True).wait()
process2=subprocess.Popen("python compile_jobs.py",shell=True).wait()
process3=subprocess.Popen("python make_submit_script.py",shell=True).wait()
