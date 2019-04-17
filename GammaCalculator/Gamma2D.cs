using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GammaCSharp
{
    class Point2D
    {
        public float X;
        public float Y;
        public Point2D()
        {

        }
        public Point2D(float x, float y)
        {
            X = x;
            Y = y;
        }
    }
    class Gamma2D
    {
        private bool global;
        private bool normalization;
        private float cutValue; //  percent of MaxDose
        private float limit;
        private float doseTol;
        private float distTol;
        private Point2D pixelSpacing = new Point2D();
        private int voxelsX;
        private int voxelsY;
        private int TotalVoxelsNum;
        private float[] refDose;
        private float[] tarDose;
        private float[] gammaIndex_;
        public float[] gammaIndex
        {
            get { return gammaIndex_; }
        }
        private string refFile;
        private string tarFile;
        private List<KeyValuePair<float, Point2D>> box = new List<KeyValuePair<float, Point2D>>();

        public Gamma2D(string refFile,string tarFile,float pixelX,float pixelY, float doseTol = 3, float distTol = 3, float cutValue = (float)0.1, bool global = true, bool normalization = false, int limit = 2)
        {
            this.refFile = refFile;
            this.tarFile = tarFile;
            pixelSpacing.X = pixelX;
            pixelSpacing.Y = pixelY;

            this.global = global;
            this.normalization = normalization;
            this.cutValue = cutValue; //  percent of MaxDose
            this.doseTol = doseTol;
            this.distTol = distTol;

            this.limit = limit;
            import();
            RunGammaCal();
            export();
        }
        public Gamma2D(float[] refDose, float[] tarDose,int voxelsX, int voxelsY, float pixelX, float pixelY, float doseTol = 3, float distTol = 3, float cutValue = (float)0.1, bool global = true, bool normalization = false, int limit = 2)
        {
            this.refDose = refDose;
            this.tarDose = tarDose;

            float refDoseMax = this.refDose.Max();
            if (normalization==true)
            {
                for (int i = 0; i < refDose.Count(); i++)
                {
                    this.refDose[i] = this.refDose[i] / refDoseMax;
                    this.tarDose[i] = this.tarDose[i] / refDoseMax;
                }
            }

            pixelSpacing.X = pixelX;
            pixelSpacing.Y = pixelY;

            this.global = global;
            this.normalization = normalization;
            this.cutValue = cutValue; //  percent of MaxDose
            this.doseTol = doseTol;
            this.distTol = distTol;
            this.limit = limit;
            this.voxelsX =voxelsX;// Numbers of Voxels X
            this.voxelsY =voxelsY;
            TotalVoxelsNum = voxelsX * voxelsY;
            refDose = new float[TotalVoxelsNum];
            tarDose = new float[TotalVoxelsNum];
            gammaIndex_ = new float[TotalVoxelsNum];

            RunGammaCal();
            export();
        }
        private void import()
        {
            string[] refdose = File.ReadAllLines(refFile);
            string[] tardose = File.ReadAllLines(tarFile);
            voxelsY = refdose.Count();
            voxelsX = refdose[0].TrimStart().TrimEnd().Split().Where(item => item != "").ToArray().Count();

            TotalVoxelsNum = voxelsX * voxelsY;
            refDose = new float[TotalVoxelsNum];
            tarDose = new float[TotalVoxelsNum];
            gammaIndex_ = new float[TotalVoxelsNum];
            int flag = 0;
            for (int i = 0; i < voxelsY; i++)
            {
                string[] refdoseline = refdose[i].TrimStart().TrimEnd().Split().Where(item => item != "").ToArray();
                string[] tardoseline = tardose[i].TrimStart().TrimEnd().Split().Where(item => item != "").ToArray();
                for (int j = 0; j < voxelsX; j++)
                {
                    refDose[flag] = Convert.ToSingle(refdoseline[j]);
                    tarDose[flag] = Convert.ToSingle(tardoseline[j]);
                    flag++;
                }
            }
        }
        private void RunGammaCal()
        {
            float maxDose = refDose.Max();

            CreateSearchBox(distTol*limit,pixelSpacing.X/ 50);
            
            for (int i = 0; i < TotalVoxelsNum; i++)
            {
                gammaIndex_[i] = limit;
                if (refDose[i] < maxDose * cutValue)
                {
                    gammaIndex_[i] = -1;
                    continue;
                }
                int[] voxelIndex = index2sub(i);
                for (int j = 0; j < box.Count; j++)
                {
                    Point2D pos = new Point2D();
                    pos.X = voxelIndex[0] * pixelSpacing.X + box[j].Value.X;
                    pos.Y = voxelIndex[1] * pixelSpacing.Y + box[j].Value.Y;
                    if (pos.X < 0 || pos.Y<0|| pos.X > (voxelsX - 1) * pixelSpacing.X|| pos.Y > (voxelsY - 1) * pixelSpacing.Y)
                        continue;

                    float refdose = My2DInterpolation(refDose, pos);
                    float distCrittol = 0;
                    if (global==true)
                        distCrittol = doseTol / 100 * maxDose;
                    else
                        distCrittol = doseTol / 100 * refdose;


                    float tardose = tarDose[i];

                    float dosediff = refdose - tardose;

                    float gammaSq = GammaSquare(dosediff * dosediff, box[j].Key, distTol * distTol, distCrittol * distCrittol);
                    float gamma =(float)Math.Sqrt(gammaSq);
                    gammaIndex_[i] =Math.Min(gamma, gammaIndex_[i]);

                    //if (box[j].Key/(distTol* distTol) > gammaIndex[i])
                    //    break;
                }

            }
            for (int i = 0; i < TotalVoxelsNum; i++)
            {
                Console.WriteLine("Gamma[" + i + "]:" + gammaIndex_[i]);
            }
        }
        private void export()
        {
            int dayu1 = gammaIndex_.ToList().Where(item => item > 1).ToList().Count();
            int dayu0 = gammaIndex_.ToList().Where(item => item > 0).ToList().Count();

            Console.WriteLine("大于0:" + dayu0);
            Console.WriteLine("大于1:" + dayu1);
            Console.WriteLine("GammaPassRate:" + (dayu0 - dayu1) * 100.0 / dayu0 + "%");
        }
        private float GammaSquare(float doseDiffSq, float distSq, float distTolSq, float doseTolSq)
        {
            return doseDiffSq / doseTolSq + distSq / distTolSq;
        }
        private float My2DInterpolation(float[] dose, Point2D point1)
        {
            Point2D point = new Point2D(point1.X / pixelSpacing.X, 
                                        point1.Y / pixelSpacing.Y);

            int X_floor = (int)Math.Floor(point.X);
            int Y_floor = (int)Math.Floor(point.Y);
            float alpha = point.X - X_floor;
            float beta = point.Y - Y_floor;
            if (alpha == 0&& beta==0)
                return dose[sub2index(X_floor,Y_floor)];
            int px, py;
            px = py = 1;
            if (Math.Abs(point.X - (voxelsX - 1)) < 1e-5) px = 0;
            if (Math.Abs(point.Y - (voxelsY - 1)) < 1e-5) py = 0;

            float XL_YL = dose[sub2index(X_floor,Y_floor)];
            float XL_YR = dose[sub2index(X_floor,Y_floor+py)];
            float XR_YL = dose[sub2index(X_floor+px,Y_floor)];
            float XR_YR = dose[sub2index(X_floor+px,Y_floor+py)];

            float value =    alpha * beta * XR_YR
                            + alpha * (1 - beta) * XR_YL
                            + (1 - alpha) * beta * XL_YR
                            + (1 - alpha) * (1 - beta) * XL_YL;
            return value;

        }
        private int[] index2sub(int index)
        {
            return new int[] {  index % voxelsX,
                                index / voxelsX % voxelsY};
        }
        private int sub2index(int x,int y)
        {
            return x + y * voxelsX;

        }
        private void CreateSearchBox(float radius, float pixel)
        {
            int voxelsNum = Convert.ToInt32(radius / pixel);
            for (int i = -voxelsNum; i <= voxelsNum; i++)
            {
                for (int j = -voxelsNum; j <= voxelsNum; j++)
                {
                    Point2D point = new Point2D();
                    point.X = i * pixel;
                    point.Y = j * pixel;
                    box.Add(new KeyValuePair<float, Point2D>(point.X * point.X + point.Y * point.Y, point));

                }
            }
            box = box.Where(item=>item.Key<=radius*radius).OrderBy(item => item.Key).ToList();
        }
    }
 }
