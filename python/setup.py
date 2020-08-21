from distutils.core import setup, Extension

ext = Extension(name='sampen', 
                sources=['sampenmodule.cpp'], 
                extra_compile_args=['-std=c++11'], 
                libraries=["sampen"], 
                library_dirs=["/home/phree/local/lib"], 
                include_dirs=['/home/phree/local/include'])

setup(name="sampen", 
      version="1.0", 
      description="Aimed at Fast Entropy Calculation. ", 
      author="phree", 
      author_email="liu.wf@outlook.com", 
      ext_modules=[ext])
