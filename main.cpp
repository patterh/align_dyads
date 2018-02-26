//============================================================================
// Name        : align_dyads.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "main.h"

using namespace std;

int main()
{
bool percentage = false;
bool log_two = true;
// Read in file of list of filenames of nucleosome dyad data that you want to analyze
// It is assumed that the order of the data files are the same as the chromosome order
// You have to load the data files for all the chromosomes, even if you don't need them all

HP_ReadTextFile* file = new HP_ReadTextFile;
char* dyad_filename = new char[1024];
strcpy(dyad_filename, "/home/hpatterton/Siegel_EMBOJ/list_of_dyad_files_Siegel.txt");
char* outfilename = new char[1024];
strcpy(outfilename, "/home/hpatterton/Siegel_EMBOJ/aligned_dyads_Siegel_log2_av.txt");
HP_DynamicStringArray* list_of_filenames = new HP_DynamicStringArray;
int number_of_files = file->ReadFileofFileNames(dyad_filename, list_of_filenames);


// Read in nucleosome dyad data files
float** array_of_dyad_values = new float*[number_of_files];
int* number_of_dyad_values = new int[number_of_files];
for(int file_number = 0; file_number < number_of_files; file_number++)
{
number_of_dyad_values[file_number] = file->ReadNucleosomeDyadFile(array_of_dyad_values[file_number], list_of_filenames->GetStringPointer(file_number));
cout << "Number of dyad values: " << number_of_dyad_values[file_number] << endl;
}
float total_dyad_value = 0;
float average_dyad_value;
int non_zero_dyad_values = 0;
int total_number_of_dyad_positions = 0;
for(int file_number = 0; file_number < number_of_files; file_number++)
{
	total_number_of_dyad_positions += number_of_dyad_values[file_number];
	for(int positions = 0; positions < number_of_dyad_values[file_number]; positions++)
	{
		total_dyad_value += array_of_dyad_values[file_number][positions];
		if(array_of_dyad_values[file_number][positions])
			non_zero_dyad_values+=1;
	}
}
average_dyad_value = total_dyad_value/total_number_of_dyad_positions;
cout << "Total dyad value: " << total_dyad_value << endl;
cout << "Total number of dyad positions: " << total_number_of_dyad_positions << endl;
cout << "Average dyad value: " << average_dyad_value << endl;
cout << "Non_zero_dyad_values: " << non_zero_dyad_values << endl;

//Calculate standard deviation

float sum_of_square_differences = 0;
float sum_of_square_differences_normalised = 0;
float sum_of_square_differences_normalised_squareroot = 0;
float difference;
float square_of_difference;
for(int file_number = 0; file_number < number_of_files; file_number++)
{
	for(int positions = 0; positions < number_of_dyad_values[file_number]; positions++)
	{
		difference = array_of_dyad_values[file_number][positions] - average_dyad_value;
		square_of_difference = difference*difference;
		sum_of_square_differences += square_of_difference;
	}
}
sum_of_square_differences_normalised = sum_of_square_differences/(total_number_of_dyad_positions-1);
sum_of_square_differences_normalised_squareroot = sqrt(sum_of_square_differences_normalised);

float* temp3 = new float[number_of_dyad_values[8]];
for(int x = 0; x < number_of_dyad_values[8]; x++)
	temp3[x] = log2(array_of_dyad_values[8][x]/average_dyad_value);
for(int x = 0; x <100; x++)
	cout << temp3[x] << endl;

float standard_deviation = file->StandardDeviation(temp3, number_of_dyad_values[8]);
cout << "Standard deviation: " << standard_deviation << endl;



// Read in file of chromosome number and genomic positions to align

char* genomic_data_filename = new char[1024];
strcpy(genomic_data_filename, "/home/hpatterton/Siegel_EMBOJ/SAS divergent.txt");
HPDynamicIntArray* chromosome_number = new HPDynamicIntArray;
HPDynamicIntArray* genomic_position = new HPDynamicIntArray;
HP_DynamicStringArray* strand = new HP_DynamicStringArray;
int number_of_genomic_features = file->ReadFileofGenomicPositions(chromosome_number, genomic_position, strand, genomic_data_filename);
cout << "features = " << number_of_genomic_features << endl;

// Retrieve relevant regions from read-in data files and copy to float array

int upstream_range = 500;
int downstream_range = 500;
int section_size = upstream_range + downstream_range + 1; // the +1 is due to the nucleotide in the center of the two ranges
float** dyad_section = new float*[number_of_genomic_features];
for(int x = 0; x < number_of_genomic_features; x++)
	{
	dyad_section[x] = new float[section_size];
	memset(dyad_section[x], 0, section_size*sizeof(float));
	}
int start;
int end;
for(int x = 0; x < number_of_genomic_features; x ++)
	{
	cout << x << ": " << genomic_position->GetEntry(x) << endl;
	if(strand->GetStringPointer(x)[0] == 'W')
		{
		//int hp = chromosome_number->GetEntry(x)-1;
		start = genomic_position->GetEntry(x)-upstream_range;
		end = genomic_position->GetEntry(x)+downstream_range;
		if((start < 1) | (end > number_of_dyad_values[chromosome_number->GetEntry(x)-1]))
			cout << "Query feature " << x << " is out of range" << endl;
		else
			{
			for(int y = 0; y < section_size; y++)
				{
				dyad_section[x][y] = array_of_dyad_values[(chromosome_number->GetEntry(x)-1)][start+y-1];
				}
			}
		}
	else
		{
		start = genomic_position->GetEntry(x)+upstream_range;
		end = genomic_position->GetEntry(x)-downstream_range;
		//int hp = chromosome_number->GetEntry(x)-1;
		if((end < 1) | (start > number_of_dyad_values[(chromosome_number->GetEntry(x)-1)]))
			cout << "Query feature " << x << " is out of range" << endl;
		else
			{
			for(int y = 0; y < section_size; y++)
				{
				dyad_section[x][y] = array_of_dyad_values[(chromosome_number->GetEntry(x)-1)][start-y-1];
				}
			}
		}
	}


// Calculate sum of float arrays

float* average_dyad_values = new float[section_size];
float sum;
for(int x = 0; x < section_size; x++)
	{
	sum = 0;
	for(int y = 0; y < number_of_genomic_features; y++)
		{
		if(dyad_section[y][x] >= 0) // test to see if we wrote a '-1' in the data file at this position, and ignore entry if true
			{
			sum += dyad_section[y][x];
			}
		}
	average_dyad_values[x] = sum;
	}

//convert to log2 if selected
float* ratio_dyad_values = new float[section_size];
if(log_two)
{
for(int x = 0; x < section_size; x++)
	{
	ratio_dyad_values[x] = log2((average_dyad_values[x]/number_of_genomic_features)/average_dyad_value);
	}
}

//Perform a running average

int running_average_range = 50;
for(int x = 0; x <= section_size-running_average_range; x++)
	{
	sum = 0.0;
	for(int y = x; y < x+running_average_range; y++)
		sum += ratio_dyad_values[y];
	ratio_dyad_values[x] = sum/running_average_range;
	}

// See if we want to express number of dyads as percentage of total
if(percentage)
{
	float total = 0;
	for(int y = 0; y <= section_size-running_average_range; y++)
		total += average_dyad_values[y];
	for(int y = 0; y <= section_size-running_average_range; y++)
		average_dyad_values[y] = 100*average_dyad_values[y]/total;
	cout << total << endl;
}



// Write result file to disk


fstream outfile;
outfile.open(outfilename, fstream::out | fstream::binary);
char* temp = new char[1024];
for(int x = 0; x < section_size; x++)
{
	memset(temp, 0, 1024);
	sprintf(temp, "%.5f%c%c", ratio_dyad_values[x],13,10);
	outfile << temp;
}
outfile.close();
cout << section_size << " values written to " << outfilename << endl;
cout << "done" << endl;

// Reclaim memory

delete [] temp;
delete [] dyad_filename;
delete [] genomic_data_filename;
for(int x = 0; x < number_of_files; x++)
	delete [] array_of_dyad_values[x];
for(int x = 0; x < number_of_genomic_features; x++)
	delete [] dyad_section[x];
delete [] dyad_section;
delete [] array_of_dyad_values;
delete [] number_of_dyad_values;
delete [] average_dyad_values;

return 0;
}
