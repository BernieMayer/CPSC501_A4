#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <math.h>

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE   16

/*  Standard sample size in bytes  */
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)

struct AudioFileHeader readHeaderOfAudioFile(FILE* file);
void testWrite();
void convolve(float x[], int N, float h[], int M, float y[], int P);
long getFileSize(FILE* file);
void readFile(float fileData[], int size, FILE* file, int offset);
void readFileDataIntoArray(short inputData[],int sizeOfInputData,  struct AudioFileHeader header, FILE* file);
void normalizeArray(float array[], int size);
void convertShortArrayToFloat(short short_Array[], int size,  float float_Array[]);

void divideArrayByItsCurrentMax(float array[], int size);
void denormalizeArray(float array[], int size);
void writeToWaveFile(float data[], int channels, int numberSamples, double outputRate, FILE* file);

size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);


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
		array[i] = array[i]/(SHRT_MAX);
	}

}

//This function will divide the array by the current absolute maximum of the
//array
void divideArrayByItsCurrentMax(float array[], int size)
{
  int max = array[0];   //Assume user will have a non zero sized array
  for (int i = 0; i < size; i ++)
  {
    if (abs(array[i]) > max)
    {
      max = abs(array[i]);
    }
  }

  for (int i = 0; i < size; i++)
  {
    array[i] = array[i]/max;
  }
}

//This fucntion will multiply each element in the float
// array by 65535 which is the maximum
void denormalizeArray(float array[], int size)
{
  for (int i = 0; i < size; i++)
  {
    array[i] = array[i] * (SHRT_MAX -1); // SHRT_MAX -1 is for safe measure so that the numbers don't overflow
  }
}

//This function will write to the wav file
// the header and the data
void writeToWaveFile(float data[], int channels, int numberSamples, double outputRate, FILE* file)
{
  writeWaveFileHeader(channels, numberSamples, outputRate, file); //use the writeWaveFileHeader function to write to the file
  //write to the file the data content
}






/******************************************************************************
*
*       function:       writeWaveFileHeader
*
*       purpose:        Writes the header in WAVE format to the output file.
*
*       arguments:      channels:  the number of sound output channels
*                       numberSamples:  the number of sound samples
*                       outputRate:  the sample rate
*                       outputFile:  the output file stream to write to
*
*       internal
*       functions:      fwriteIntLSB, fwriteShortLSB
*
*       library
*       functions:      ceil, fputs
*
******************************************************************************/

void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;

    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;

    /*  Calculate the total number of bytes per frame  */
    short int frameSize = channels * BYTES_PER_SAMPLE;

    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);

    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);

    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);

    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((short)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}



/******************************************************************************
*
*       function:       fwriteIntLSB
*
*       purpose:        Writes a 4-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}

/******************************************************************************
*
*       function:       fwriteShortLSB
*
*       purpose:        Writes a 2-byte integer to the file stream, starting
*                       with the least significant byte (i.e. writes the int
*                       in little-endian form).  This routine will work on both
*                       big-endian and little-endian architectures.
*
*       internal
*       functions:      none
*
*       library
*       functions:      fwrite
*
******************************************************************************/

size_t fwriteShortLSB(short int data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
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
