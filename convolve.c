#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>

struct AudioFileHeader readHeaderOfAudioFile(FILE* file);
void testWrite();
void convolve(float x[], int N, float h[], int M, float y[], int P);
long getFileSize(FILE* file);
void readFile(float fileData[], int size, FILE* file, int offset);
void readFileDataIntoArray(float inputData[],int sizeOfInputData,  struct AudioFileHeader header, FILE* file);
void normalizeArray(float array[], int size);

struct  AudioFileHeader
{
  int chunkId;
  int chunkSize;
  int format;

  int subchunk1ID;
  int subchunk1Size;
  unsigned short audioFormat;
  unsigned short numChannels;
  int sampleRate;
  int byteRate;
  unsigned short blockAlign;
  unsigned short bitsPerSample;
  int subchunk2ID;
  unsigned int subchunk2Size;
};


void readFile(float fileData[], int size, FILE* file, int offset)
{
  fseek(file, offset, SEEK_SET);

  int i;

  
}

struct AudioFileHeader readHeaderOfAudioFile(FILE* file)
{
  int j = 0;
  int i = 0;
  struct AudioFileHeader audioHeader;
  //fscanf(file, "%d", &i);

  fread(&audioHeader, sizeof(struct AudioFileHeader), 1, file);
  return audioHeader;
}


//This convolve function is based on the one in the CPSC 501 lecture slides
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
  int n,m;

  /* Clear output buffer y[] */
  for ( n = 0; n < P; n++)
    y[n] = 0.0;
	
    /*Outer loop process each input value x[n] */
    for ( n = 0; n <N; n++) {
	  //printf("The convolution is %f%% done \n", ((float) n/ (float) N) * 100);
      //Inner loop process x[n] with each sample of h[n]
      for (m = 0; m <M; m++)
      {
        y[n+m] += x[n] + h[m];
      }
    }
}

// gets the file size in bytes
long getFileSize(FILE* file)
{
  //fopen(file);
  fseek(file, 0, SEEK_END);
  long result = ftell(file);
  rewind(file);
  //fclose(file)
  printf("The file size has been calculated");
  return result;
}

void testWrite()
{
  int write = 134;

  printf("About to make the test file");

  FILE* file = fopen("tmp","w+" );
  printf("Test file has been created");

  fwrite(&write, sizeof(write), 1, file );
  fclose(file);

  printf("Test file has been made");

}

void readFileDataIntoArray(float inputData[], int sizeOfInputData, struct AudioFileHeader header, FILE* file)
{
	int chunkSize = header.subchunk2Size;
	int bitsPerSample = header.bitsPerSample;
	int frameSize = header.blockAlign;
	
	int i;
	
	fseek(file, sizeof(header), SEEK_SET);
	printf("Just about to read the file data into the array \n");
	for (i = 0; i < sizeOfInputData; i++)
	{
		fread(&inputData[i], bitsPerSample, 1, file);
	}
	
	
	
	
}

void normalizeArray(float array[], int size)
{
	int i;
	
	for (i = 0; i < size; i++)
	{
		array[i] = array[i]/FLT_MAX;
	}
	
}
	
int main(int argc, char * argv[])
{


  /*
  Testing splitting an int
  int num = 0b00110010;
  unsigned short part1 = num & 0xF;
  unsigned short part2 = (num >> 4) &0xF;
  printf("The num is %d , part1 is %hu part2 is %hu", num, part1, part2);
  */
  
  long startTime = time(NULL);


  FILE* inputFile;
  FILE* IRFile;
  FILE* outputFile;
  if (argc > 2)
  {
	  //This block of code initializes the inputFile data
      char* inputFileName = argv[1];
      printf("The file name is %s \n" , inputFileName);
      inputFile = fopen(inputFileName, "r");	//Open the input file
      
      //The AudioFileHeader takes data from the header of the file and uses this data 
      //as the program see fit
      // this data is stored in a structure.
      struct AudioFileHeader inputAudioHeader = readHeaderOfAudioFile(inputFile);
	  int sizeOfInputData = inputAudioHeader.subchunk2Size/(inputAudioHeader.bitsPerSample * inputAudioHeader.numChannels); // Also is the number of Samples
      float inputData[sizeOfInputData];
      readFileDataIntoArray(inputData, sizeOfInputData, inputAudioHeader, inputFile);
     
      
      //This block of code initializes the IR file data
      
      char* IRFileName = argv[2];
      IRFile = fopen(IRFileName, "r");
      
      struct AudioFileHeader IR_AudioHeader = readHeaderOfAudioFile(IRFile);
      //printf("Done making the IR header \n");
      int numSamplesIR = IR_AudioHeader.subchunk2Size/(IR_AudioHeader.bitsPerSample * IR_AudioHeader.numChannels);
	  int sizeOfIRData = 1000000;
      
      float IR_Data[sizeOfIRData];
      
      //readFileDataIntoArray(inputData, sizeOfInputData, inputAudioHeader, inputFile);
      readFileDataIntoArray(IR_Data, sizeOfIRData, IR_AudioHeader, IRFile);
      
      //printf("done reading the IR file data \n");
      
      
      //normalize the arrays to be between 1 and -1..
      normalizeArray(inputData, sizeOfInputData);
      normalizeArray(IR_Data, sizeOfIRData);
      
      int sizeOutput = sizeOfInputData + sizeOfIRData -1;
      float outputData[sizeOutput];
      printf("Starting the convolution \n");
      
      
      int sizeLeft = numSamplesIR - sizeOfIRData;
      while (sizeLeft > 0) {
		printf("The convolution is %f%% done \n", ( (float) sizeLeft / (float) numSamplesIR  * 100));
		convolve(inputData, sizeOfInputData, IR_Data, sizeOfIRData, outputData, sizeOutput); 
		//will need to convolve with all of the IR_Data... 
		sizeLeft = sizeLeft - sizeOfIRData;
      }
      
      if (sizeLeft < 0) {
		  convolve(inputData, sizeOfInputData, IR_Data, numSamplesIR, outputData, sizeOutput); 
	  }
      printf("Done the convolution \n");
      
      long finishTime = time(NULL);
      
      printf("The convolution took %li ", finishTime - startTime);
      
      
      //Handle the output file creation here...
      

  }


  return 0;
}
