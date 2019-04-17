using Dicom;
using Dicom.Imaging;
using Dicom.Imaging.Mathematics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GammaCSharp
{
    class Gamma3D
    {

        private Stopwatch stopwatch = new Stopwatch();

        private Point3D pixelSpacing = new Point3D();
        private int voxelsX; // no. of voxels in x-direction
        private int voxelsY;
        private int voxelsZ;
        private float doseTol; // [Gy]
        private float distTol; // [mm]
        private bool global;   // default global == 1 to calculate global gamma index
        private bool normalization;
        private int cutoffCount;
        private float cutValue;
        private float limit;

        private int totalVoxelNum; 
        private float[] refDose;
        private float[] tarDose;
        public  float[] gammaIndex;
        private string refRTDoseFile;
        private string tarRTDoseFile;

        private List<KeyValuePair<float, Point3D>> box = new List<KeyValuePair<float, Point3D>>();


        public Gamma3D(string refRTDoseFile, string tarRTDoseFile,  float doseTol = 3, float distTol = 3, float cutValue = (float)0.1, bool global = true, bool normalization = false, int limit = 2)
        {
            this.refRTDoseFile = refRTDoseFile;
            this.tarRTDoseFile = tarRTDoseFile;
            this.global = global;
            this.normalization = normalization;
            this.cutValue = cutValue; //  percent of MaxDose
            this.doseTol = doseTol; //in percent
            this.distTol = distTol; //[mm]
            this.limit = limit;

            Import();
            stopwatch.Start();
            RunGammaCal();
            stopwatch.Stop();
            Export();
        }

        public Gamma3D(float[] refDose, float[] tarDose,int[] Sizes,float[] pixelSpacings, float doseTol = 3, float distTol = 3, float cutValue = (float)0.1, bool global = true, bool normalization = false, int limit = 2)
        {
            this.refDose = refDose;
            this.tarDose = tarDose;
            // dose Normalization
            float refDoseMax = this.refDose.Max();
            if (normalization == true)
            {
                for (int i = 0; i < totalVoxelNum; i++)
                {
                    refDose[i] = this.refDose[i]  / refDoseMax;
                    tarDose[i] = this.tarDose[i]  / refDoseMax;
                }
            }

            this.global = global;

            this.normalization = normalization;
            this.cutValue = cutValue; //  percent of MaxDose
            this.doseTol = doseTol; //in percent
            this.distTol = distTol; //[mm]
            this.limit = limit;

            voxelsX = Sizes[0];
            voxelsY = Sizes[1];
            voxelsZ = Sizes[2];
            totalVoxelNum = voxelsX * voxelsY * voxelsZ;
            refDose = new float[totalVoxelNum];
            tarDose = new float[totalVoxelNum];
            gammaIndex = new float[totalVoxelNum];
            pixelSpacing.X = pixelSpacings[0];
            pixelSpacing.Y = pixelSpacings[1];
            pixelSpacing.Z = pixelSpacings[2];

            stopwatch.Start();
            RunGammaCal();
            stopwatch.Stop();
            Export();
        }

        private void Import()
        {
            DicomDataset dt1 = DicomFile.Open(refRTDoseFile).Dataset;
            DicomDataset dt2 = DicomFile.Open(tarRTDoseFile).Dataset;

            float[] PixelSpacing = dt1.GetValues<float>(DicomTag.PixelSpacing);
            float SliceThickness = dt1.GetSingleValue<float>(DicomTag.SliceThickness);
            voxelsX = dt1.GetSingleValue<int>(DicomTag.Columns);
            voxelsY = dt1.GetSingleValue<int>(DicomTag.Rows);
            voxelsZ = dt1.GetSingleValue<int>(DicomTag.NumberOfFrames);
            totalVoxelNum = voxelsX * voxelsY * voxelsZ;
            refDose = new float[totalVoxelNum];
            tarDose = new float[totalVoxelNum];
            gammaIndex = new float[totalVoxelNum];

            pixelSpacing.X = PixelSpacing[0];
            pixelSpacing.Y = PixelSpacing[1];
            pixelSpacing.Z = SliceThickness;
            UInt16[] dt1dosedata = dt1.GetValues<UInt16>(DicomTag.PixelData);
            float DoeGridScaling1 = dt1.GetSingleValue<float>(DicomTag.DoseGridScaling);
            UInt16[] dt2dosedata = dt2.GetValues<UInt16>(DicomTag.PixelData);
            float DoeGridScaling2 = dt2.GetSingleValue<float>(DicomTag.DoseGridScaling);

            float maxvalue1 = dt1dosedata.Max() * DoeGridScaling1;

            // dose Normalization
            if(normalization == true)
            {
                for (int i = 0; i < totalVoxelNum; i++)
                {
                    refDose[i] = (dt1dosedata[i] * DoeGridScaling1) / maxvalue1;
                    tarDose[i] = (dt2dosedata[i] * DoeGridScaling2) / maxvalue1;
                }
            }
            else
            {
                for (int i = 0; i < totalVoxelNum; i++)
                {
                    refDose[i] = (dt1dosedata[i] * DoeGridScaling1);
                    tarDose[i] = (dt2dosedata[i] * DoeGridScaling2);
                }
            }

            //TextWriter tw = new StreamWriter("d://text.txt");
            //for (int i = 0; i < voxelsY; i++)
            //{
            //    for (int j = 0; j < voxelsX; j++)
            //    {
            //        int x = Sub2Index(j, i, 100);
            //        tw.WriteLine(refDose[x] + ",Sub:(" + j + "," + i + ",100) Index:" + x);
            //        //tw.Write(refDose[x]+" ");
            //        float dd = 1;
            //    }
            //    //tw.WriteLine();
            //}
            //tw.Flush();


        }

        /// <summary>
        /// 对Eva进行重采样
        /// </summary>
        /// <param name="interpretNum">点数扩大倍数</param>
        private void ImportWithResample(int interpretNum)
        {
            DicomDataset dt1 = DicomFile.Open(refRTDoseFile).Dataset;
            DicomDataset dt2 = DicomFile.Open(tarRTDoseFile).Dataset;

            // Ref读取
            float[] PixelSpacing = dt1.GetValues<float>(DicomTag.PixelSpacing);
            float SliceThickness = dt1.GetSingleValue<float>(DicomTag.SliceThickness);
            voxelsX = dt1.GetSingleValue<int>(DicomTag.Columns);
            voxelsY = dt1.GetSingleValue<int>(DicomTag.Rows);
            voxelsZ = dt1.GetSingleValue<int>(DicomTag.NumberOfFrames);
            totalVoxelNum = voxelsX * voxelsY * voxelsZ;
            refDose = new float[totalVoxelNum];
            gammaIndex = new float[totalVoxelNum];

            // Eva读取
            tarDose = new float[totalVoxelNum];
            pixelSpacing.X = PixelSpacing[0];
            pixelSpacing.Y = PixelSpacing[1];
            pixelSpacing.Z = SliceThickness;
            UInt16[] dt1dosedata = dt1.GetValues<UInt16>(DicomTag.PixelData);
            float DoeGridScaling1 = dt1.GetSingleValue<float>(DicomTag.DoseGridScaling);
            UInt16[] dt2dosedata = dt2.GetValues<UInt16>(DicomTag.PixelData);
            float DoeGridScaling2 = dt2.GetSingleValue<float>(DicomTag.DoseGridScaling);

            float maxvalue1 = dt1dosedata.Max() * DoeGridScaling1;

            // dose Normalization
            if (normalization == true)
            {
                for (int i = 0; i < totalVoxelNum; i++)
                {
                    refDose[i] = (dt1dosedata[i] * DoeGridScaling1) / maxvalue1;
                    tarDose[i] = (dt2dosedata[i] * DoeGridScaling2) / maxvalue1;
                }
            }
            else
            {
                for (int i = 0; i < totalVoxelNum; i++)
                {
                    refDose[i] = (dt1dosedata[i] * DoeGridScaling1);
                    tarDose[i] = (dt2dosedata[i] * DoeGridScaling2);
                }
            }

            //TextWriter tw = new StreamWriter("d://text.txt");
            //for (int i = 0; i < voxelsY; i++)
            //{
            //    for (int j = 0; j < voxelsX; j++)
            //    {
            //        int x = Sub2Index(j, i, 100);
            //        tw.WriteLine(refDose[x] + ",Sub:(" + j + "," + i + ",100) Index:" + x);
            //        //tw.Write(refDose[x]+" ");
            //        float dd = 1;
            //    }
            //    //tw.WriteLine();
            //}
            //tw.Flush();


        }

        private void RunGammaCal()
        {
            float maxDose = refDose.Max();
            Console.WriteLine("maxDose = " + maxDose);
            CreateSearchBox(distTol*limit, (float)pixelSpacing.X/20);
            cutoffCount = refDose.ToList().Where(item=>item<maxDose*cutValue).ToList().Count();

            //Parallel.For(0, totalVoxelNum, (idx, loopState1) =>
            for(int i=0;i<totalVoxelNum;i++)
            {
                gammaIndex[i] = limit;
                float DistCritSq ;
                float DoseCritSq ;

                if (refDose[i] < maxDose * cutValue)
                {
                    gammaIndex[i] = -1;
                    continue;
                }
                //Determine x,y,z from idx of the reference vol
                int[] voxelIndex = Index2Sub(i);
                for (int j = 0; j < box.Count(); j++)
                {
                    Point3D SearchPos = new Point3D(
                                         voxelIndex[0] * pixelSpacing.X + box[j].Value.X,
                                         voxelIndex[1] * pixelSpacing.Y + box[j].Value.Y,
                                         voxelIndex[2] * pixelSpacing.Z + box[j].Value.Z);
                    float SearchPosRefDose;
                    if ((SearchPos.X < 0) || (SearchPos.Y < 0) || (SearchPos.Z < 0) || 
                                     (SearchPos.X > (voxelsX * pixelSpacing.X - 1)) || 
                                     (SearchPos.Y > (voxelsY * pixelSpacing.Y - 1)) || 
                                     (SearchPos.Z > (voxelsZ * pixelSpacing.Z - 1)))
                    {
                        SearchPosRefDose = 0;
                    }
                    else
                        SearchPosRefDose = MyInterpolate(refDose, SearchPos);

                    //float SearchPosTarDose = MyInterpolate(tarDose, SearchPos);
                    float tardose = tarDose[i];

                    DistCritSq = box[j].Key / (distTol * distTol);

                    ////dose evaluated within searchbox
                    if (global == true)
                        DoseCritSq = (float)Math.Pow((tardose - SearchPosRefDose) / (maxDose * doseTol / 100.0), 2);//global
                    else
                        DoseCritSq = (float)Math.Pow((tardose - SearchPosRefDose) / (SearchPosRefDose * doseTol / 100.0), 2);//local

                    float gammaSq = DistCritSq + DoseCritSq;
                    if (gammaSq < gammaIndex[i])
                    {
                        gammaIndex[i] = gammaSq;// gamma^2
                        if (gammaSq <= 1)
                            break;
                    }
                    //if the distance is larger than that the gamma for that point then break
                    if (DistCritSq > gammaIndex[i])
                        break;
                }
                if (i % 1000 == 1)
                {
                    Console.WriteLine("gammaIndex[ " + i + "] = " + gammaIndex[i]);
                }
            }//);
        }

        private void Export()
        {
            int dayu1 = gammaIndex.Where(item => item >= 1.0).ToList().Count();
            int dayu0 = gammaIndex.Where(item => item >= 0).ToList().Count();

            Console.WriteLine("花费了" + stopwatch.Elapsed.TotalSeconds + "秒");
            Console.WriteLine("总像素数:" + totalVoxelNum);
            Console.WriteLine("CUTOFF:" + cutoffCount);
            Console.WriteLine("计算个数:" + (totalVoxelNum- cutoffCount).ToString());
            Console.WriteLine("大于0个数："+dayu0);
            Console.WriteLine("大于1个数:" + dayu1);
            Console.WriteLine("gammaPassRate:"+((dayu0-dayu1) / Convert.ToSingle(dayu0) *100).ToString("f4")+"%");
        }
 
        private int Sub2Index(int x, int y, int z)
        {
           return x   + y * voxelsX + z * voxelsX * voxelsY;
        }
        private int[] Index2Sub(int index)
        {
            return new int[] {  index % voxelsX,
                                index / voxelsX % voxelsY,
                                index / (voxelsX * voxelsY) % voxelsZ};
        }
        //Inverse Distance Weighted Interpolate
        private float MyInterpolate(float[] refDose, Point3D point1)
        {
            Point3D point = new Point3D( point1.X / pixelSpacing.X, 
                                         point1.Y / pixelSpacing.Y, 
                                         point1.Z / pixelSpacing.Z );
            // Cooridinate to Index
            int X_floor = (int)Math.Floor(point.X);
            int Y_floor = (int)Math.Floor(point.Y);
            int Z_floor = (int)Math.Floor(point.Z);
            // Weights
            float alpha =(float) (point.X - X_floor);
            float beta = (float)(point.Y - Y_floor);
            float gamma = (float)(point.Z - Z_floor);
            if (alpha == 0 && beta == 0 && gamma == 0)
                return refDose[Sub2Index(X_floor,Y_floor,Z_floor)];

            int px,py,pz;
            px = py = pz = 1;
            if (Math.Abs(point.X - (voxelsX - 1)) < 1e-5) px = 0;
            if (Math.Abs(point.Y - (voxelsY - 1)) < 1e-5) py = 0;
            if (Math.Abs(point.Z - (voxelsZ - 1)) < 1e-5) pz = 0;

            //get dose of the 8 peak voxels
            float XL_YL_ZL = refDose[Sub2Index(X_floor , Y_floor , Z_floor)];               
            float XL_YL_ZR = refDose[Sub2Index(X_floor , Y_floor , Z_floor + pz)];           
            float XL_YR_ZL = refDose[Sub2Index(X_floor , Y_floor + py , Z_floor )];       
            float XL_YR_ZR = refDose[Sub2Index(X_floor , Y_floor + py , Z_floor + pz)];   
            float XR_YL_ZL = refDose[Sub2Index(X_floor + px , Y_floor , Z_floor )];       
            float XR_YL_ZR = refDose[Sub2Index(X_floor + px , Y_floor , Z_floor + pz)];    
            float XR_YR_ZL = refDose[Sub2Index(X_floor + px , Y_floor + py, Z_floor )];
            float XR_YR_ZR = refDose[Sub2Index(X_floor + px , Y_floor + py, Z_floor + pz)];

            return    alpha * beta * gamma * XR_YR_ZR
                    + alpha * beta * (1 - gamma) * XR_YR_ZL
                    + alpha * (1 - beta) * gamma * XR_YL_ZR
                    + alpha * (1 - beta) * (1 - gamma) * XR_YL_ZL
                    + (1 - alpha) * beta * gamma * XL_YR_ZR
                    + (1 - alpha) * beta * (1 - gamma) * XL_YR_ZL
                    + (1 - alpha) * (1 - beta) * gamma * XL_YL_ZR
                    + (1 - alpha) * (1 - beta) * (1 - gamma) * XL_YL_ZL;

        }
        private void CreateSearchBox(float radius, float boxStep)
        {
            int boxWidth =Convert.ToInt32(radius/boxStep);

            for (int i = -boxWidth; i <= boxWidth; i++)
                for (int j = -boxWidth; j <= boxWidth; j++)
                    for (int k = -boxWidth; k <= boxWidth; k++)
                    {
                        Point3D point = new Point3D(i * boxStep, j * boxStep, k * boxStep);
                        float distanceSq =(float)( point.X * point.X + point.Y * point.Y + point.Z * point.Z);
                        box.Add(new KeyValuePair<float, Point3D>(distanceSq, point));
                    }
            float radiusSq = radius * radius;
            box = box.Where(item=>item.Key<= radiusSq).ToList().OrderBy((item => item.Key)).ToList();
            int a = 0;
            //TextWriter tw = new StreamWriter("d://text.txt");
            //foreach (var d in box)
            //{
            //    if (d.Key < radius)
            //    {
            //        tw.WriteLine(d.Value);
            //    }
            //}
        }
    }
}
