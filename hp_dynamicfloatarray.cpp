#include "hp_dynamicfloatarray.h"

HPDynamicFloatArray::HPDynamicFloatArray(int size, int number)
{
delta = number;
upperlimit = size;
array = new float[upperlimit];
members = 0;
}

HPDynamicFloatArray::~HPDynamicFloatArray()
{
delete [] array;
}

HPDynamicFloatArray::HPDynamicFloatArray(const HPDynamicFloatArray& passedobject)
{
delta = passedobject.delta;
upperlimit = passedobject.upperlimit;
members = passedobject.members;
array = new float[passedobject.upperlimit];
for(int x = 0; x < passedobject.members; x++)
	array[x] = passedobject.array[x];

}

int HPDynamicFloatArray::AddFloat(float number)
{

if(members+1 > upperlimit)
	{
	float* temp = new float[members];
	for(int x = 0; x < members; x++)
		temp[x] = array[x];
	delete [] array;
	array = new float[members+delta];
	upperlimit = members+delta;
	for(int y = 0; y < members; y++)
		array[y] = temp[y];
	array[members] = number;
	delete [] temp;
	}
else
	{
	array[members] = number;
	}
return(++members);
}

int HPDynamicFloatArray::GetArrayLength()
{
return(members);
}

float* HPDynamicFloatArray::GetArray()
{
return(array);
}

float HPDynamicFloatArray::GetEntry(int number)
{
if(number < members)
	return(array[number]);
else
	return((float)members);
}

float HPDynamicFloatArray::GetMaximum()
{
float maximum = array[0];
for(int x = 1; x < members; x++)
	maximum = (maximum > array[x] ? maximum : array[x]);
return(maximum);
}

float HPDynamicFloatArray::GetMinimum()
{
float minimum = array[0];
for(int x = 1; x < members; x++)
	minimum = (minimum < array[x] ? minimum : array[x]);
return(minimum);
}

float HPDynamicFloatArray::GetRange()
{
return(GetMaximum()-GetMinimum());
}
