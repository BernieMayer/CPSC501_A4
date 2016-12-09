#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <math.h>

//MACRO functions
#define SWAP(a,b) tempr = (a); (a) = (b); (b) = tempr

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
void normalizeArrayDouble(double array[], int size);
void convertShortArrayToFloat(short short_Array[], int size,  float float_Array[]);
void convertShortArrayToDouble(short short_Array[], int size, double double_Array[]);

void divideArrayByItsCurrentMax(float array[], int size);
void denormalizeArray(float array[], int size);
void writeToWaveFile(float data[], int channels, int numberSamples, double outputRate, FILE* file);
void writeWaveFileHeader(int channels, int numberSamples, double outputRate, FILE *outputFile);
void writeWavFileContent(float data[], int numberSamples, FILE* file );
size_t fwriteIntLSB(int data, FILE *stream);
size_t fwriteShortLSB(short int data, FILE *stream);
void print_vector(char *title, float x[], int N);
void testIdentityConvole();
void convertFloatArrayToShort(float float_Array[], int size, short short_Array[]);
void writeToWaveFile(float data[], int channels, int numberSamples, double outputRate, FILE* file);
//void writeWavFileContent(short data_short[], int channels, int numberSamples, double outputRate, FILE *file);


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
  int i = 0;

  for (int i = 0; i < size; i++)
  {
    float_Array[i] = (float) short_Array[i];
  }

}

void convertShortArrayToDouble(short short_Array[], int size, double double_Array[])
{
  int i = 0;

  for (int i = 0; i < size; i++)
  {
    double_Array[i] = (double) short_Array[i];
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

//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
    double data_1 = data[j];
    double data_2 = data[j + 1];
		tempr = wr * data_1 - wi * data_2;
		tempi = wr * data_2 + wi * data_1;
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}


//Convolve function using Fast Fourier Transform
//(might want to convoert to doubles)
void convolveFFT(double x[], int N, double h[], int M, double y[], int P)
{

  //make segments
    //Note: the segments should be close to the size of h[]

  printf("Starting convolution using FFT \n");

  int SEGMENT_SIZE = pow(2, 14);

  if (M > SEGMENT_SIZE)
  {
    printf("The segment size needs to be changed since it will not fill the entire IR \n");
  }

  int numSegmentsX = 2;

  double * padded_h = (double * ) malloc(sizeof(double) * 2 * SEGMENT_SIZE);
  double * segmentArray_x = (double * ) malloc(sizeof(double) * 2 * SEGMENT_SIZE);


  //take all of h into padded_h and for each element in h add a 0 to pad it
  for (int j = 0; j < 2 * SEGMENT_SIZE; j++)
  {
    //Pad the h array so that it is SEGMENT_SIZE
    if (j < M)
      padded_h[j] = (double) h[j];
    else
      padded_h[j] = (double) 0.0f;
  }

  for (int i = 0; i < P; i++)
    y[i] = 0.0;


    four1(padded_h - 1,SEGMENT_SIZE, 1);

  //convolute the arrays
  double * results = (double *) malloc(sizeof(double) * SEGMENT_SIZE * 2);
  int baseIndex = 0;
  for (int i = 0; i < numSegmentsX; i++)
  {
    for (int j = 0; j < 2 * SEGMENT_SIZE; j++)
    {
      if ( j < SEGMENT_SIZE)
        segmentArray_x[j] = x[baseIndex + j];
      else
        segmentArray_x[j] = 0.0;
    }
    //baseIndex += SEGMENT_SIZE - 1;
    four1(segmentArray_x - 1, SEGMENT_SIZE, 1);

    for (int k = 0; k < SEGMENT_SIZE; k+=2)
    {
      float a = segmentArray_x[k];  //real component
      float b = segmentArray_x[k + 1];  //imaginary component

      float c = padded_h[k];  //real component
      float d = padded_h[k+1];  //imaginary compenent

      results[k] = a * c - b * d;
      results[k + 1] = c * b + a * d;



      four1(results -1, SEGMENT_SIZE, -1);
    }
    for (int r = 0; r < SEGMENT_SIZE * 2; r+=2)
    {
      if ((baseIndex + r) < P && (baseIndex + 1 + r) < P) {

        y[baseIndex + r] += (float) results[r];
        y[baseIndex + r + 1] += (float) results[r + 1];
      }
    }
    baseIndex += SEGMENT_SIZE - 1;
  }
  //for (int i = 0; i < numSegmentsX; i++)
  //{
    //Apply the fast fourier transform on the array
    //note the code below assumes an array starting at 1.
    // By subtracting by 1 we are getting the right data
    //four1(segmentArray_x - 1,2* SEGMENT_SIZE, 1);

    //Apply the FFT to the IR data
    //four1(padded_h, 2* SEGMENT_SIZE, 1);

    //convolve the data...
    //Reminder (a + b * i) * (c + d * i) = a * c + c * b * i + a * d * i - b * d
    // where i is an imaginary number
	/*
    for (int j = 0; j < SEGMENT_SIZE; j+=2)
    {
      float a = segmentArray_x[i][j];
      float b = segmentArray_x[i][j +1]; //imaginary coefficient

      float c = padded_h[j];
      float d = padded_h[j + 1];

      float real_coefficientOutput = a * c - b * d;
      float imag_coefficientOutput = c * b + a * d;

      segmentArray_y[i][j] = real_coefficientOutput;
      segmentArray_y[i][j] = imag_coefficientOutput;

      //Re convolve the data back to normal
      four1(segmentArray_y[i] - 1, SEGMENT_SIZE, -1);
    }
	*/
  //}
	printf("Why ?? \n");

	//free(padded_h);

  }








//This convolve function is based on the one in the CPSC 501 lecture slides
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
  int n,m;

  printf("N is %i and M is %i and P is %i \n", N, M,P);

  /* Clear output buffer y[] */
  for ( n = 0; n < P; n++)
    y[n] = 0.0;

    /*Outer loop process each input value x[n] */
    for ( n = 0; n <N; n++) {
	    //printf("The convolution is %f%% done \n", ((float) n/ (float) N) * 100);
      //Inner loop process x[n] with each sample of h[n]
      for (m = 0; m <M; m++)
      {
        y[n+m] += x[n] * h[m];
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
		fread(&inputData[i], sizeof(short), 1, file);
	}




}

void normalizeArray(float array[], int size)
{
	int i;

	for (i = 0; i < size; i++)
	{

    float temp = (float)array[i]/( (float)(SHRT_MAX));
    if (temp > 1.0)
    {
      array[i] = 1.0;
    } else {
      array[i] = -1.0;
    }
    array[i] = temp;

	}

}

void normalizeArrayDouble(double array[], int size)
{
  int i;

  for (i = 0; i < size; i++)
  {

    double temp = (double)array[i]/( (double)(SHRT_MAX));
    if (temp > 1.0)
    {
      array[i] = 1.0;
    } else {
      array[i] = -1.0;
    }
    array[i] = temp;

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
  printf("Max is %u", max);
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
  printf("done making the header \n");
  writeWavFileContent(data, numberSamples, file);
}

void writeWavFileContent(float data[], int numberSamples,FILE* file )
{
  int maximumValue = (int)pow(2.0, (double)BITS_PER_SAMPLE - 1) - 1;
  for (int i = 0;  i < numberSamples; i++)
  {
    double value = data[i];
    short int sampleValue = rint(value * 0.5f * maximumValue);
    //printf("sampleValue is %f");
    fwriteShortLSB(sampleValue, file);
  }

}

void convertFloatArrayToShort(float float_Array[], int size, short short_Array[])
{
  for (int i = 0; i < size; i++)
  {
    short_Array[i] = (short)( float_Array[i] * (SHRT_MAX -1));
  }
}
/*

Likely will need removal for production

void writeWavFileContent(short data_short[], int channels, int numberSamples, double outputRate, FILE *file)
{/*
  for (int i = 0; i < numberSamples; i++)
  {

  }
  int numberOfSamples = numberSamples;

  float PI = 3.14159265358979;

  int maximumValue = (int)pow(2.0, (double)BITS_PER_SAMPLE - 1) - 1;
  double angularFrequency = 2.0 * PI * 440.0;
  double increment = angularFrequency / SAMPLE_RATE;
  for (int i = 0; i < numberOfSamples; i++) {
      /*  Calculate the sine wave in the range -1.0 to + 1.0  */
      //double value = data_short[i];

      /*  Convert the value to a 16-bit integer, with the
          range -maximumValue to + maximumValue.  The calculated
          value is rounded to the nearest integer  */
      //short int sampleValue = rint(value * maximumValue);

      /*  Write out the sample as a 16-bit (short) integer
          in little-endian format  */
      //fwriteShortLSB(sampleValue, file);
  //}
//}






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

void testConvolve()
{
  float input_signal[100], impulse_response[20], output_signal[120];
  int input_size, impulse_size, output_size;

  /*  Create an example input signal  */
  input_signal[0] = 1.0;
  input_signal[1] = 0.5;
  input_signal[2] = 0.25;
  input_signal[3] = 0.125;
  input_size = 4;



  /*  Create an "identity" impulse response.  The output should be
      the same as the input when convolved with this  */
  impulse_response[0] = 1.0;
  impulse_size = 1;

  output_size = input_size + impulse_size - 1;

  convolveFFT(input_signal, input_size, impulse_response, impulse_size,
     output_signal, output_size);

  int result = 1;

  for (int i = 0; i < output_size; i++)
  {
    if (output_signal[i] != input_signal[i])
      printf("Convolution test failed, expected %f but got %f", input_signal[i], output_signal[i]);
  }
  if (result)
  {
    printf("Convolve test has succeeded \n");
  }
}

void testConvolve2()
{

  float input_signal[100], impulse_response[20], output_signal[120];
  int input_size, impulse_size, output_size;

  /*  Create an example input signal  */
  input_signal[0] = 1.0;
  input_signal[1] = 0.5;
  input_signal[2] = 0.25;
  input_signal[3] = 0.125;
  input_size = 4;

  /*  Create an "echo effect".  The output will contain the original signal
      plus a copy delayed by 2 samples and 1/2 the amplitude.  The original
      and copy will overlap starting at the 3rd sample  */
  impulse_response[0] = 1.0;
  impulse_response[1] = 0.0;
  impulse_response[2] = 0.5;
  impulse_size = 3;

  /*  Set the expected size of the output signal  */
  output_size = input_size + impulse_size - 1;

  /*  Do the convolution, and print the output signal  */
  convolve(input_signal, input_size, impulse_response, impulse_size,
     output_signal, output_size);



    print_vector("echo",output_signal, output_size);
}


void print_vector(char *title, float x[], int N)
{
  int i;

  printf("\n%s\n", title);
  printf("Vector size:  %-d\n", N);
  printf("Sample Number \tSample Value\n");
  for (i = 0; i < N; i++)
    printf("%-d\t\t%f\n", i, x[i]);
}

void testProfile()
{
  printf("Profiling the program  3.0 \n");
  FILE* inputFile;
  FILE* IRFile;
  FILE* outputFile;


  inputFile = fopen("GuitarDry.wav", "r");
  IRFile = fopen("l960taj_mahal_deep6.wav", "r");
  outputFile = fopen("profile.wav", "w+");

  //run the base program ...
  for (int i = 0;  i < 10; i++) {
    printf("Running the program for the %i th time", i);
    //seek the files to the be used
    fseek(inputFile, 0, SEEK_SET);
    fseek(IRFile, 0, SEEK_SET);
    fseek(outputFile, 0, SEEK_SET);

    struct AudioFileHeader inputAudioHeader = readHeaderOfAudioFile(inputFile);

    unsigned int sizeOfInputData = inputAudioHeader.subchunk2Size/( (float)(inputAudioHeader.bitsPerSample/ 8) * inputAudioHeader.numChannels); // Also is the number of Samples
    short * inputDataShort  = (short *) malloc(sizeof(short) * sizeOfInputData);


    //printf("The number of bits per sample is %i", inputAudioHeader.bitsPerSample);
    //printf("The size of a short is %i", sizeof(short));
    //printf("The size of subchunk2Size is %i", inputAudioHeader.subchunk2Size);


    readFileDataIntoArray(inputDataShort, sizeOfInputData, inputAudioHeader, inputFile);

    double* inputData = (double *) malloc( sizeof(double) * sizeOfInputData);
    convertShortArrayToDouble(inputDataShort, sizeOfInputData, inputData);
    normalizeArray(inputData, sizeOfInputData);

    struct AudioFileHeader ir_AudioHeader = readHeaderOfAudioFile(IRFile);

    unsigned int sizeOfIRData = ir_AudioHeader.subchunk2Size/(2.0f);  //Should be non magic number
    short* ir_Data_short = (short *) malloc(sizeof(short) * sizeOfIRData);
    double* ir_data = (double *) malloc(sizeof(double) * sizeOfIRData);
    convertShortArrayToDouble(ir_Data_short, sizeOfIRData, ir_data);


    int sizeOutput = sizeOfInputData + sizeOfIRData - 1;
    double * outputData = (double *) malloc(sizeof(double) * sizeOutput);

    convolveFFT(inputData, sizeOfInputData, ir_data, sizeOfIRData, outputData, sizeOutput);
    int test = 1;
    for (int i = 0; i < 10; i++)
    {
      if (inputData[i] != outputData[i])
      {
        printf("test failed expected %f but got %f... on i is %i \n", inputData[i], outputData[i], i);
        test = 0;
      }
    }

    if (test == 1)
    {
      printf("Test has succeeded \n");
    }

    outputFile = fopen("testID.wav", "w");
    writeToWaveFile(inputData, 1, sizeOutput, (double) inputAudioHeader.sampleRate, outputFile);


    //free up any memory used


    free(ir_Data_short);
    free(ir_data);
  }

  fclose(inputFile);
  fclose(IRFile);
  fclose(outputFile);

}



int main(int argc, char * argv[])
{
    testProfile();

    return 0;
    //testIdentityConvole();




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
	  int sizeOfInputData = inputAudioHeader.subchunk2Size/((float)(inputAudioHeader.bitsPerSample/8) * inputAudioHeader.numChannels); // Also is the number of Samples
      //int sizeOfInputData = 100;

      if (debug) {
      printf("The sizeOfInputData is %i \n", sizeOfInputData);
      printf("The size of subchunk1D is  %u", inputAudioHeader.subchunk1Size);
      printf("The number of bits per sample is %u \n", inputAudioHeader.bitsPerSample);
      printf("The byteRate is  %u \n", inputAudioHeader.byteRate);
      printf("The subchunk2ID is %u \n", inputAudioHeader.subchunk2ID);
      printf("The size of the input file is %u bytes \n", inputAudioHeader.subchunk2Size);
      }
      short* inputDataShort = (short *) malloc(sizeof(short) *sizeOfInputData);
      //float inputData[sizeOfInputData];
      readFileDataIntoArray(inputDataShort, sizeOfInputData, inputAudioHeader, inputFile);



      double* inputData = (double * ) malloc(sizeof(double) * sizeOfInputData);

      convertShortArrayToDouble(inputDataShort, sizeOfInputData, inputData);

      //This block of code initializes the IR file data

      char* IRFileName = argv[2];
      printf("The IRFile nams is %s \n", IRFileName);
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
      double IR_Data[numSamplesIR];




      readFileDataIntoArray(IR_Data_short, numSamplesIR, inputAudioHeader, inputFile);


      convertShortArrayToDouble(IR_Data_short, numSamplesIR, IR_Data);



      printf("done reading the IR file data \n");





      //normalize the arrays to be between 1 and -1..
      normalizeArrayDouble(inputData, sizeOfInputData);
      normalizeArrayDouble(IR_Data, numSamplesIR);

      if (debug) {
      for (int i = 0; i < sizeOfInputData; i++)
      {
        if (inputData[i] > 1.0 || inputData[i] < -1.0)
        {
          inputData[i] = 1.0;
          printf("normalization failed... inputData[%u] is %f \n", i ,inputData[i] );

        }
      }
    }


      int sizeOutput = sizeOfInputData + numSamplesIR -1;
      printf("The sizeOutput is %u", sizeOutput);
      double* outputData = (double *) malloc(sizeof(double) *sizeOutput);



      printf("Starting the convolution \n");
      printf("sizeOutput is %u while sizeOfInputData is %u and numSamplesIR is %u", sizeOutput, sizeOfInputData, numSamplesIR);
      convolveFFT(inputData, sizeOfInputData, IR_Data, numSamplesIR, outputData, sizeOutput);

      long finishTime = time(NULL);



      printf("The convolution took %li ", finishTime - startTime);

      char* outputFileName = argv[3];
      outputFile = fopen( outputFileName, "w");

      printf("The output file name is %s \n" , outputFileName);



      divideArrayByItsCurrentMax(outputData, sizeOutput);
      //denormalizeArray(outputData, sizeOutput);
      printf("The data has been renormalized");
      if (debug)
      {
        //print_vector("Data has been renormalized",outputData, sizeOutput);
      }

      //print_vector("Just a test...", outputData, sizeOutput);

      writeToWaveFile(outputData, 1, sizeOutput, (double) IR_AudioHeader.sampleRate, outputFile);





      printf("Done the convolution \n");

      //fclose(inputFile);
      //fclose(IRFile);
      //fclose(outputFile);



      free(inputDataShort);
      free(outputData);

      //Handle the output file creation here...


  }


  return 0;
}
