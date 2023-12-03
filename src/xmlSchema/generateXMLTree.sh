#remove old files
rm schema.cpp
rm schema.h
#generate the new ones
xsdcxx cxx-tree --std c++11 schema.xsd
#rename the new files
mv schema.cxx schema.cpp
mv schema.hxx schema.h
#change include line from <#include "schema.hpp"> to <#include "schema.h">
sed -i '41s/.*/#include "schema.h"/' schema.cpp