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
    Array<DEM::Particle *> p;        // the array of particles at which the force is to be applied;
    Array<DEM::Particle *> e;        // the array of particles which are fixed;
    String             test;         // Type of test vibraiton or tension
    double             A;            // Area of the plate for stress calculation
    double             Am;           // vibration amplitude
    double             ome;          // vibration frequency
    double             sy;           // Stress State
    double             Tf;           // Final time
    double             ex;           // Finlat strain
    double             theta;        // Angle of deflection
    Vec3_t             L0;           // Initial dimensions
    std::ofstream      oss_ss;       // file for stress strain data
};

void Setup (DEM::Domain & Dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dat.p.Size()>1)
    {
        double ome    = dat.ome*4.0*M_PI/dat.Tf*(tanh((Dom.Time-0.5*dat.Tf)/(0.1*dat.Tf)));
        for (size_t ip=0;ip<dat.p.Size();ip++)
        {
            dat.p[ip]-> w (2)  = 0.5*ome;
            dat.p[ip]-> wb(2)  = 0.5*ome;
            dat.e[ip]-> w (2)  = 0.5*ome;
            dat.e[ip]-> wb(2)  = 0.5*ome;
        }
        if (Dom.Time>0.5*dat.Tf)
        {
            for (size_t ip=0;ip<dat.p.Size();ip++)
            {
                dat.p[ip]->vxf = false;
                dat.p[ip]->vyf = false;
                dat.p[ip]->vzf = false;
                dat.p[ip]->wxf = false;
                dat.p[ip]->wyf = false;
                //dat.p[ip]->wzf = false;

                dat.e[ip]->vxf = false;
                dat.e[ip]->vyf = false;
                dat.e[ip]->vzf = false;
                dat.e[ip]->wxf = false;
                dat.e[ip]->wyf = false;
                //dat.e[ip]->wzf = false;
            }
            double Am = dat.Am*tanh((Dom.Time-0.5*dat.Tf)/(0.1*dat.Tf));
            double r0 = Dom.Particles[0]->Dmax*2.0/sqrt(3.0);
            #pragma omp parallel for schedule(static) num_threads(Dom.Nproc)
            for (size_t n=0;n<Dom.ListPosPairs.Size();n++)
            {
                size_t i = Dom.ListPosPairs[n].first;
                size_t j = Dom.ListPosPairs[n].second;
                bool pi_has_vf = !Dom.Particles[i]->IsFree();
                bool pj_has_vf = !Dom.Particles[j]->IsFree();
                bool same_tag  = Dom.Particles[i]->Tag==Dom.Particles[j]->Tag;
                double r = DEM::Distance(Dom.Particles[i]->x,Dom.Particles[j]->x);
                if ( pi_has_vf || pj_has_vf || (r>2*r0) || (r<r0) || same_tag) continue;
                double r6 = pow((r0/r),6.0);
                Vec3_t n = (Dom.Particles[i]->x - Dom.Particles[j]->x)/r;
                Vec3_t F = Am/r*r6*(1.0-r6)*n;

                omp_set_lock  (&Dom.Particles[i]->lck);
                Dom.Particles[i]->F -= F;
                omp_unset_lock(&Dom.Particles[i]->lck);
                omp_set_lock  (&Dom.Particles[j]->lck);
                Dom.Particles[j]->F += F;
                omp_unset_lock(&Dom.Particles[j]->lck);
            }
        }
    }
    else 
    {
        double ome = dat.ome*6.0*M_PI/dat.Tf*(tanh((Dom.Time-dat.Tf/3.0)/(0.1*dat.Tf)));
        if (Dom.Time>dat.Tf/3.0) ome = 0.0;
        dat.p[0]->w (2)  = ome;
        dat.p[0]->wb(2)  = ome;
        if (Dom.Time>dat.Tf/3.0&&Dom.Time<=2.0*dat.Tf/3.0)
        {
            dat.p[0]->FixVeloc(3.0*dat.ex*dat.L0(1)/dat.Tf*sin(dat.theta*M_PI/180.0),3.0*dat.ex*dat.L0(1)/dat.Tf*cos(dat.theta*M_PI/180.0),0.0);
            dat.p[0]->InitializeVelocity(Dom.Dt);
        }
        else if (Dom.Time>2.0*dat.Tf/3.0)
        {
            dat.p[0]->FixVeloc(0.0,0.0,0.0);
            dat.p[0]->InitializeVelocity(Dom.Dt);
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
    size_t RenderVideo; // Decide if video should be rendered
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
    size_t Ns;          // Number of strings
    double Am;          // vibration force amplitude or angle of deflection
    double ome;         // Frequency of vibration or number of torsion cycles
    double ex;          // Final strain for the tensile test (positive extension, negative compression)
    double theta;       // Angle of deflection for strain tests
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
        infile >> Ns;           infile.ignore(200,'\n');
        infile >> Am;           infile.ignore(200,'\n');
        infile >> ome;          infile.ignore(200,'\n');
        infile >> ex;           infile.ignore(200,'\n');
        infile >> theta;        infile.ignore(200,'\n');
    }



    // user data and domain
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha = verlet;
    //dom.Dilate= true;

    if (ptype=="voronoi") 
    {
        for(size_t ns=0;ns<Ns;ns++)
        {
            size_t Initialidx = dom.Particles.Size();
            dom.AddVoroPack (-3*ns, R, Lx,Ly,Lz, nx,ny,nz, rho, true, false, seed, 1.0, Vec3_t(0.9,0.9,0.9));
            for (size_t ip=Initialidx;ip<dom.Particles.Size();ip++)
            {
                dom.Particles[ip]->Translate(Vec3_t(0.0,0.0,ns*sqrt(2.0)*Lz));
            }
            dom.GenBoundingPlane(-3*ns-1,-3*ns,R,1.0,true);
            dat.p.Push(dom.GetParticle(-3*ns-2));
            dat.e.Push(dom.GetParticle(-3*ns-1));
            dom.GetParticle(-3*ns-1)->FixVeloc();
            dom.GetParticle(-3*ns-2)->FixVeloc();
        }
    }
    else throw new Fatal("Packing for particle type not implemented yet");

    // Initialize the UserData structure
    dat.test = test;
    dat.A    = Lx*Ly;
    dat.Am   = Am;
    dat.ome  = ome;
    dat.Tf   = Tf;
    dat.ex   = ex;
    dat.theta= theta;
    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    dat.L0   = Xmax - Xmin;


    //set the element properties
    Dict B;
    for (size_t ns=0;ns<Ns;ns++)
    {
        B.Set(-3*ns  ,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
        B.Set(-3*ns-1,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
        B.Set(-3*ns-2,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    }
    dom.SetProps(B);

    //dom.WriteXDMF("rubber");
    

    dom.Solve (Tf,dt,dtOut, &Setup, &Report, filekey.CStr(), RenderVideo, Nproc);
}
MECHSYS_CATCH
