using Dicom;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NPlot;
using System.Drawing;
using System.IO;

namespace GammaCSharp
{
    class Program_BP
    {
        static void Main_BP(string[] args)
        {
            string dir = @"D:\Sync\STUDY\Team\TOPAS\AutoTOPAS\AutoTOPAS\TOPASFile\PBS_work\gammaTest\gamma-matlab\";

            int test =2;
            if(test==1)
            {
                string ref1Dfile = dir + "ref1D.txt";
                string tar1Dfile = dir + "tar1D.txt";
                Gamma1D gamma1 = new Gamma1D(ref1Dfile, tar1Dfile, pixelSpacing: 2, cutValue:(float)0.1,doseTol: 2, distTol: 2, global: true, normalization: false, limit: 2);
            }
            else if(test==2)
            {

                string ref2Dfile = dir + "ref2D.txt";
                string tar2Dfile = dir + "tar2D.txt";
                Gamma2D gamma2 = new Gamma2D(ref2Dfile, tar2Dfile, pixelX: 2,pixelY:2, doseTol: 2, distTol: 2, cutValue:(float) 0.1, global: true, normalization: false, limit: 2);
            }
            else
            {
                string ref3Dfile = dir + "MergeBeams1RD1.2.752.243.1.1.20181009165012577.2700.51523.dcm";
                string tar3Dfile = dir + "MergeBeams2RD1.2.752.243.1.1.20181009165012577.2700.51523.dcm";
                Gamma3D gamma3 = new Gamma3D(ref3Dfile, tar3Dfile, doseTol: 3, distTol: 3, cutValue: (float)0.1, global: true, normalization: false, limit: 2);
            }

            Console.ReadKey();
        }
    }
}
