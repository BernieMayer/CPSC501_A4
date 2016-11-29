#include <stdio.h>
#include <string.h>
#include <stdlib.h>

struct AudioFileHeader readHeaderOfAudioFile(FILE* file);
void testWrite();

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


struct AudioFileHeader readHeaderOfAudioFile(FILE* file)
{
  int j = 0;
  int i = 0;
  struct AudioFileHeader audioHeader;
  //fscanf(file, "%d", &i);

  fread(&i, sizeof(int), 1, file);
  while (!feof(file))
  {
    printf("%d\n", i);
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
      audioHeader.blockAlign = i & 0xF;
      audioHeader.bitsPerSample = (i >> 4) & 0xF;
    } else if ( j == 9)
      audioHeader.subchunk2ID = i;
    else if (j == 10)
      audioHeader.subchunk2Size = i;
    fread(&i, sizeof(int), 1, file);
    j++;
  }
  fclose(file);
  return audioHeader;
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
  }


  return 0;
}
