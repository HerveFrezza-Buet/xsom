# xsom

A C++ generic programming library for multi-som experiments, by <a href="https://github.com/jeremyfix">Jérémy Fix</a> and <a href="https://github.com/HerveFrezza-Buet">Hervé Frezza-Buet</a>.

# Unix Installation

First, get the files.

``` 
git clone https://github.com/HerveFrezza-Buet/xsom
``` 

Then, you can install the package as follows. 

```
mkdir -p xsom/build
cd xsom/build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr
sudo make install
cd ../..
```

For a Fedora-64bit architecture:

```
mkdir -p xsom/build
cd xsom/build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr -DLIB_SUFFIX=64
sudo make install
cd ../..
```

