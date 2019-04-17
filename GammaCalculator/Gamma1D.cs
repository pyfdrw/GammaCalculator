using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Dicom;
using Dicom.Imaging.Mathematics;

namespace GammaCSharp
{
    class Gamma1D
    {
        private bool global;
        private bool normalization;
        private float cutValue; //  percent of MaxDose
        private float doseTol;
        private float distTol;
        private string refFile;
        private string tarFile;
        private float pixelSpacing;
        private float[] refDose;
        private float[] tarDose;
        private float[] gammaIndex;
        private int TotalVoxelsNum;
        private float limit;
        private List<KeyValuePair<float, float>> box = new List<KeyValuePair<float, float>>();

        public Gamma1D(string refFile,string tarFile,float pixelSpacing, float doseTol = 3,float distTol = 3,float cutValue = (float)0.1,bool global=true,bool normalization=false,int limit=2)
        {
            this.global = global;
            this.normalization = normalization;
            this.cutValue = cutValue; //  percent of MaxDose
            this.doseTol = doseTol;
            this.distTol = distTol;
            this.limit = limit;
            this.refFile = refFile;
            this.tarFile = tarFile;
            this.pixelSpacing = pixelSpacing;

            import();
            RunGammaCal();
            export();
        }
        private void import()
        {
            string[] refdose = File.ReadAllLines(refFile);
            string[] tardose = File.ReadAllLines(tarFile);
            TotalVoxelsNum = refdose.Count();
            refDose = new float[TotalVoxelsNum];
            tarDose = new float[TotalVoxelsNum];
            gammaIndex = new float[TotalVoxelsNum];
            for (int i = 0; i < TotalVoxelsNum; i++)
            {
                refDose[i] =Convert.ToSingle(refdose[i]);
                tarDose[i] = Convert.ToSingle(tardose[i]);
            }

            if (normalization==true)
            {
                float refDoseMax = refDose.Max();
                float tarDoseMax = tarDose.Max();
                for (int i = 0; i < TotalVoxelsNum; i++)
                {
                    refDose[i] = refDose[i] / refDoseMax;
                    tarDose[i] = tarDose[i] / tarDoseMax;

                }
            }
        }
        private void RunGammaCal()
        {
            float maxDose = refDose.Max();
            CreateSearchBox(distTol*limit, pixelSpacing/100);
            for (int i = 0; i < TotalVoxelsNum; i++)
            {
                gammaIndex[i] = limit;
                if (refDose[i] < maxDose * cutValue)
                {
                    gammaIndex[i] = -1;
                    continue;
                }
                for (int j = 0; j < box.Count; j++)
                {
                    float posX = i * pixelSpacing + box[j].Value;
                    float refdose;
                    if (posX < 0 || posX > (TotalVoxelsNum - 1) * pixelSpacing)
                        refdose = 0;
                    else
                        refdose = My1DInterpolation(refDose,posX);
                    //float tardose = My1DInterpolation(tarDose,posX);
                    float doseCritTol;
                    if (global==true)
                    {
                        doseCritTol = doseTol / 100 * maxDose;
                    }
                    else
                    {
                        doseCritTol = doseTol / 100 * refdose;
                    }

                    float tardose = tarDose[i];

                    float dosediff = refdose - tardose;

                    float gammaSq = GammaSquare(dosediff * dosediff, box[j].Key, distTol*distTol, doseCritTol * doseCritTol);
                    float gamma =(float) Math.Sqrt(gammaSq);
                    gammaIndex[i] = Math.Min(gammaIndex[i], gamma);
                }

            }
            for (int i = 0; i < TotalVoxelsNum ; i++)
            {

                Console.WriteLine("Gamma["+i+"]:" +gammaIndex[i]);
            }
        }
        private void export()
        {
            int dayu1 = gammaIndex.ToList().Where(item => item > 1).ToList().Count();
            int dayu0 = gammaIndex.ToList().Where(item => item >= 0).ToList().Count();

            Console.WriteLine("大于0:" + dayu0 );
            Console.WriteLine("大于1:" + dayu1 );
            Console.WriteLine("GammaPassRate:" + (dayu0 - dayu1) * 100.0 / dayu0+"%");
        }
        private float GammaSquare(float doseDiffSq,float distSq,float distTolSq,float doseTolSq)
        {
            return doseDiffSq / doseTolSq + distSq / distTolSq;
        }
        private float My1DInterpolation(float[] dose,float X1)
        {
            float X = X1 / pixelSpacing;
            int X_floor = (int)Math.Floor(X);
            float alpha = X - X_floor;
            if (alpha==0)
                return dose[X_floor];

            float value=(1 - alpha) * dose[X_floor] + alpha * dose[X_floor + 1];
            return value;
            
        }
        private void CreateSearchBox(float radius, float pixel)
        {
            int voxelsNum = Convert.ToInt32(radius / pixel);
            for (int i = -voxelsNum; i <= voxelsNum; i++)
            {
                float x = i * pixel;
                box.Add(new KeyValuePair<float, float>(x * x, x));
            }
            box=box.OrderBy(item=>item.Key).ToList();
        }

    }
}
