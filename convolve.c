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
void readFileDataIntoArray(short inputData[],int sizeOfInputData,  struct AudioFileHeader header, FILE* file);
void normalizeArray(float array[], int size);
void convertShortArrayToFloat(short short_Array[], int size,  float float_Array[]);

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

//This method will take in the short array and will convert it into the
// float array
void convertShortArrayToFloat(short short_Array[], int size, float float_Array[])
{

  for (int i = 0; i < size; i++)
  {
    float_Array[i] = (float) short_Array[i];
  }

}

struct AudioFileHeader readHeaderOfAudioFile(FILE* file)
{
  int j = 0;
  int i = 0;
  struct AudioFileHeader audioHeader;
  //fscanf(file, "%d", &i);

  fread(&audioHeader, sizeof(struct AudioFileHeader), 1, file);

  if (audioHeader.subchunk1Size == 18)
  {
    //This means that subchunk1 is 18 bytes..

    int subchunk2ID;
    int subchunk2Size;

    fseek(file,  (sizeof(audioHeader) - 2 * sizeof(int) + 2), SEEK_SET);
    fread(&subchunk2ID, sizeof(int), 1, file);
    fread(&subchunk2Size, sizeof(int), 1, file);

    audioHeader.subchunk2ID = subchunk2ID;
    audioHeader.subchunk2Size = subchunk2Size;

  }

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

void readFileDataIntoArray(short inputData[], int sizeOfInputData, struct AudioFileHeader header, FILE* file)
{
	int chunkSize = header.subchunk2Size;
	int bitsPerSample = header.bitsPerSample;
	int frameSize = header.blockAlign;

	int i;

	fseek(file, sizeof(header), SEEK_SET);
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
		array[i] = array[i]/(65535);
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
  int debug = 1;
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
      //int sizeOfInputData = 100;

      if (debug) {
      printf("The sizeOfInputData is %i \n", sizeOfInputData);
      printf("The size of subchunk1D is  %u", inputAudioHeader.subchunk1Size);
      printf("The number of bits per sample is %u \n", inputAudioHeader.bitsPerSample);
      printf("The byteRate is  %u \n", inputAudioHeader.byteRate);
      printf("The subchunk2ID is %u \n", inputAudioHeader.subchunk2ID);
      printf("The size of the input file is %u bytes \n", inputAudioHeader.subchunk2Size);
      }
      short inputDataShort[sizeOfInputData];
      //float inputData[sizeOfInputData];
      readFileDataIntoArray(inputDataShort, sizeOfInputData, inputAudioHeader, inputFile);

      float inputData[sizeOfInputData];

      convertShortArrayToFloat(inputDataShort, sizeOfInputData, inputData);

      //This block of code initializes the IR file data

      char* IRFileName = argv[2];
      IRFile = fopen(IRFileName, "r");

      struct AudioFileHeader IR_AudioHeader = readHeaderOfAudioFile(IRFile);
      //printf("Done making the IR header \n");
      int numSamplesIR = IR_AudioHeader.subchunk2Size/(IR_AudioHeader.bitsPerSample * IR_AudioHeader.numChannels);
	    //int sizeOfIRData =

      if (debug) {
      printf("Showing the info about the IR file \n");
      printf("The size of subchunk1D is  %u", IR_AudioHeader.subchunk1Size);
      printf("The size of the data is %u\n", IR_AudioHeader.subchunk2Size);
      printf("The byteRate is  %u \n", IR_AudioHeader.byteRate);
      printf("The subchunk2ID is %u \n", IR_AudioHeader.subchunk2ID);
      printf("The number of sameples is %u\n", numSamplesIR);
      }
      short IR_Data_short[numSamplesIR];
      float IR_Data[numSamplesIR];




      readFileDataIntoArray(IR_Data_short, numSamplesIR, inputAudioHeader, inputFile);


      convertShortArrayToFloat(IR_Data_short, numSamplesIR, IR_Data);



      printf("done reading the IR file data \n");





      //normalize the arrays to be between 1 and -1..
      normalizeArray(inputData, sizeOfInputData);
      normalizeArray(IR_Data, numSamplesIR);

      for (int i = 0; i < sizeOfInputData; i++)
      {
        if (inputData[i] > 1.0 || inputData[i] < -1.0)
        {
          printf("normalization failed...");
        }
      }




      int sizeOutput = sizeOfInputData + numSamplesIR -1;
      printf("The sizeOutput is %u", sizeOutput);
      float outputData[sizeOutput];

      convolve(inputData, sizeOfInputData, IR_Data, numSamplesIR, outputData, sizeOutput);









      printf("Done the convolution \n");

      long finishTime = time(NULL);

      printf("The convolution took %li ", finishTime - startTime);


      //Handle the output file creation here...


  }


  return 0;
}
