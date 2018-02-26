
#ifndef _HPDYNAMICFLOATARRAY_H
#define _HPDYNAMICFLOATARRAY_H

using namespace std;

class HPDynamicFloatArray
{
public:
HPDynamicFloatArray(int size = 10, int delta = 10);// size is initial size
HPDynamicFloatArray(const HPDynamicFloatArray&);     // delta is by how much array grows
~HPDynamicFloatArray();

float* GetArray();
int GetArrayLength();
float GetEntry(int index);
int AddFloat(float number);
float GetMaximum();
float GetMinimum();
float GetRange();

protected:
float* array;
int members;
int delta;
int upperlimit;

private:

HPDynamicFloatArray& operator=(const HPDynamicFloatArray&);
};

#endif
