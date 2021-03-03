import subprocess
commands='''
mkdir femdg_tutorial
cd femdg_tutorial

mkdir tmp
cd tmp
git init
git remote add origin https://gitlab.dune-project.org/dune-fem/dune-fem-dg.git
git config core.sparsecheckout true
echo 'pydemo/' > .git/info/sparse-checkout
echo 'doc/'  >> .git/info/sparse-checkout
git pull origin master
cp pydemo/camc-paper/*.py pydemo/camc-paper/*.hh pydemo/camc-paper/*.dgf ..
cd ..
rm -rf tmp
'''
subprocess.check_output(commands, shell=True)

print("###################################################################")
print("## The tutorial is now located in the 'femdg_tutorial' folder. ")
try:
    import matplotlib
except ImportError:
    print("## Note: some of the examples require the installation of 'matplotlib'.")
try:
    import scipy
except ImportError:
    print("## Note: some of the examples require the installation of 'scipy'.")
print("###################################################################")
