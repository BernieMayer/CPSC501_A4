#include <stdio.h>
#include <string.h>
#include <stdlib.h>

struct AudioFileHeader readHeaderOfAudioFile(FILE* file);
void testWrite();
void convolve(float x[], int N, float h[], int M, float y[], int P);
long getFileSize(FILE* file);
void readFile(float fileData[], int size, FILE* file, int offset);

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
  int subchunk2Size;
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
  /*
  while (!feof(file))
  {
    //printf("%d\n", i);
    if (j == 0)
      audioHeader.chunkId = i;
    else if (j == 1)
      audioHeader.chunkSize = i;
    else if (j == 2)
      audioHeader.format = i;
    else if (j == 3)
      audioHeader.subchunk1ID = i;
    else if (j == 4)
      audioHeader.subchunk1Size = i;
    else if (j == 5)
    {
      audioHeader.audioFormat = i & 0xF;
      audioHeader.numChannels = (i >> 4) & 0xF;
    } else if ( j == 6)
      audioHeader.sampleRate = i;
    else if (j == 7)
      audioHeader.byteRate = i;
    else if (j == 8)
    {
      audioHeader.blockAlign =  (short) i & 0xF;
      audioHeader.bitsPerSample = (short) (i >> 4) & 0xF;
      printf(" j is 8 i is %i \n" ,i);
    } else if ( j == 9)
      audioHeader.subchunk2ID = i;
    else if (j == 10)
      audioHeader.subchunk2Size = i;
    fread(&i, sizeof(int), 1, file);
    j++;
  }
  */
  //fclose(file);
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


int main(int argc, char * argv[])
{


  /*
  Testing splitting an int
  int num = 0b00110010;
  unsigned short part1 = num & 0xF;
  unsigned short part2 = (num >> 4) &0xF;
  printf("The num is %d , part1 is %hu part2 is %hu", num, part1, part2);
  */


  FILE* audioFile;
  if (argc > 1)
  {
      char* fileName = argv[1];
      printf("The file name is %s \n" , fileName);
      audioFile = fopen(fileName, "r");
      struct AudioFileHeader audioHeader = readHeaderOfAudioFile(audioFile);
      printf("The size of audioHeader is %i \n", sizeof(audioHeader));
      //seek past the header...

      
      printf("The size of the data chunk is %i \n", audioHeader.subchunk2Size);
      printf("The number of samples per data chunk is %hi \n", audioHeader.bitsPerSample);
      printf("The number of channels is %hi", audioHeader.numChannels);

      float inputData[audioHeader.subchunk2Size];

      //readFileDataIntoArray(inputData, size, audioFile, sizeof(audioFile));

  }


  return 0;
}
