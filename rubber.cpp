/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// Std lib
#include <math.h>

// GSL
#include <gsl/gsl_linalg.h>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;

struct UserData
{
    DEM::Particle *    p;            // the array of particles at which the force is to be applied
    Array<Vec3_t  >    vm0;          // value of the vectors close to the middle section
    Array<Vec3_t *>    vm;           // pointers to the vectors close to the middle section
    String             test;         // Type of test vibraiton or tension
    double             A;            // Area of the plate for stress calculation
    double             Am;           // vibration amplitude
    double             ome;          // vibration frequency
    double             sy;           // Stress State
    double             Tf;           // Final time
    double             ex;
    Vec3_t             L0;           // Initial dimensions
    std::ofstream      oss_ss;       // file for stress strain data
};

void Setup (DEM::Domain & Dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dat.test=="torsion")
    {
        if (Dom.Time>dat.Tf/3.0&&Dom.Time<=2.0*dat.Tf/3.0)
        {
            dat.p->w(2)   = 0.0;
            dat.p->wb(2)  = 0.0;
            //dat.p->vyf = true;
            //dat.p->v(1) = 2.0*dat.ex*dat.L0(1)/dat.Tf;
            //dat.p->xb(1) = dat.p->x(1)-dat.p->v(1)*Dom.Dt;
            dat.p->FixVeloc(3.0*dat.ex*dat.L0(1)/dat.Tf*sin(dat.Am*M_PI/180.0),3.0*dat.ex*dat.L0(1)/dat.Tf*cos(dat.Am*M_PI/180.0),0.0);
            dat.p->InitializeVelocity(Dom.Dt);
            //dat.p->wxf  = false;
            //dat.p->wyf  = false;
        }
        else if (Dom.Time>2.0*dat.Tf/3.0)
        {
            dat.p->FixVeloc(0.0,0.0,0.0);
            dat.p->InitializeVelocity(Dom.Dt);
        }
        
    }
}

void Report (DEM::Domain & Dom, void * UD)
{ 
    UserData & dat = (*static_cast<UserData *>(UD));
}

int main(int argc, char **argv) try
{

    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    size_t Nproc = 1; 
    if (argc==3) Nproc=atoi(argv[2]);
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());

    double verlet;      // Verlet distance for optimization
    String ptype;       // Particle type 
    String test;       // Particle type 
    size_t RenderVideo; // Decide is video should be render
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Bn;          // Cohesion normal stiffness
    double Bt;          // Cohesion tangential stiffness
    double Bm;          // Cohesion torque stiffness
    double Eps;         // Threshold for breking bonds
    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double dt;          // Time step
    double dtOut;       // Time step for output
    double Tf;          // Final time for the test
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    double rho;         // rho
    double Am;          // vibration force amplitude
    double ome;         // Frequency of vibration
    double ex;          // Final strain for the tensile test (positive extension, negative compression)
    {
        infile >> verlet;       infile.ignore(200,'\n');
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> test;         infile.ignore(200,'\n');
        infile >> RenderVideo;  infile.ignore(200,'\n');
        infile >> Kn;           infile.ignore(200,'\n');
        infile >> Kt;           infile.ignore(200,'\n');
        infile >> Gn;           infile.ignore(200,'\n');
        infile >> Gt;           infile.ignore(200,'\n');
        infile >> Mu;           infile.ignore(200,'\n');
        infile >> Bn;           infile.ignore(200,'\n');
        infile >> Bt;           infile.ignore(200,'\n');
        infile >> Bm;           infile.ignore(200,'\n');
        infile >> Eps;          infile.ignore(200,'\n');
        infile >> R;            infile.ignore(200,'\n');
        infile >> seed;         infile.ignore(200,'\n');
        infile >> dt;           infile.ignore(200,'\n');
        infile >> dtOut;        infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
        infile >> Lx;           infile.ignore(200,'\n');
        infile >> Ly;           infile.ignore(200,'\n');
        infile >> Lz;           infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> rho;          infile.ignore(200,'\n');
        infile >> Am;           infile.ignore(200,'\n');
        infile >> ome;          infile.ignore(200,'\n');
        infile >> ex;           infile.ignore(200,'\n');
    }



    // user data and domain
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha = verlet;
    //dom.Dilate= true;

    if (ptype=="voronoi") dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, true, false, seed, 1.0, Vec3_t(0.9,0.9,0.9));
    else if (ptype=="cube")
    {
        Mesh::Structured mesh(3);
        mesh.GenBox (false, nx, ny, nz, Lx, Ly, Lz);
        dom.GenFromMesh (mesh, R, rho, true, false);
    }
    else if (ptype=="tetra")
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
        dom.GenFromMesh (mesh, R, rho, true, false);
    }
    else throw new Fatal("Packing for particle type not implemented yet");

    // Initialize the UserData structure
    dat.test = test;
    dat.A    = Lx*Ly;
    dat.Am   = Am;
    dat.ome  = ome;
    dat.Tf   = Tf;
    dat.ex   = ex;
    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    dat.L0   = Xmax - Xmin;

    dom.GenBoundingPlane(-2,R,1.0,true);

    //identify the moving lid
    dat.p = dom.GetParticle (-3);
    if (test=="tensile")
    {
        dat.p->FixVeloc();
        dat.p->v = 0.0, ex*dat.L0(1)/Tf, 0.0;
    }


    //set the element properties
    Dict B;
    B.Set(-1,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(-2,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(-3,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    dom.SetProps(B);

    // fix -2 particles at the left extreme of the beam
    DEM::Particle * p;
    p = dom.GetParticle (-2);
    p->FixVeloc();

    if (test=="torsion")
    {
        dat.p->FixVeloc();
        //dat.p->wxf  = false;
        //dat.p->wyf  = false;
        //dat.p->wzf  = true;
        dat.p->w(2) = 6*M_PI*ome/Tf;
        //p->wxf = false;
        //p->wyf = false;
    }

    //dom.WriteXDMF("rubber");
    //

    dom.Solve (Tf,dt,dtOut, &Setup, &Report, filekey.CStr(), RenderVideo, Nproc);
}
MECHSYS_CATCH
