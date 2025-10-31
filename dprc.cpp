#ifdef DPRC

#include <cstdio>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <deepmd/deepmd.hpp>
#include "f90stream.hpp"

//
// Amber uses -DMPI when compiling with MPI, but this can't be
// defined within CXX code because it clashes with the MPI
// namespace within mpi.h
//
// Undefine MPI if it was defined, and define WITH_MPI instead
//
#ifdef MPI
#ifndef WITH_MPI
#define WITH_MPI
#endif
#undef MPI
#endif

#ifdef WITH_MPI
#include <mpi.h>
#endif




#define BOXL 200.
namespace ml
{
   static double QMCUT = 10.;

#ifdef WITH_MPI
  MPI_Comm comm;
  int mpirank=0;
  int mpisize=0;
#endif
}


void get_avg_and_std
( bool calcavg,
  std::vector<double> const & all_ene,
  std::vector< std::vector<double> > const & all_frc,
  double & ene,
  std::vector<double> & frc,
  double & maxstd )
{
  int nmodel = all_ene.size();
  ene = 0.;
  for ( int i=0; i<nmodel; ++i )
    ene += all_ene[i];
  ene /= nmodel;

  int ncrd = all_frc[0].size();
  frc.resize(ncrd);
  for ( int i=0; i<ncrd; ++i )
    frc[i] = 0.;

  int nat = ncrd / 3;

  for ( int i=0; i<ncrd; ++i )
    {
      for ( int k=0; k<nmodel; ++k )
	frc[i] += all_frc[k][i];
      frc[i] /= nmodel;
    }

  maxstd = 0.;
  for ( int a=0; a<nat; ++a )
    {
      double stddev = 0.;
      for ( int i=0; i<nmodel; ++i )
	{
	  double dx[3] = { all_frc[i][0+a*3] - frc[0+a*3],
			   all_frc[i][1+a*3] - frc[1+a*3],
			   all_frc[i][2+a*3] - frc[2+a*3] };
	  //std::printf("%4i %3i %15.6e\n",a,i,dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	  stddev += dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]; 
	}
      //std::printf("%4i %15.6e\n",a,stddev);
      stddev = std::sqrt( stddev / nmodel );
      if ( stddev > maxstd )
	maxstd = stddev;
    }

  if ( not calcavg )
    {
      ene = all_ene[0];
      for ( int i=0; i<ncrd; ++i )
	frc[i] = all_frc[0][i];
    }
}


std::vector<int> MakeIdxMap(int nat)
{
  std::vector<int> idxmap(nat);
  for ( int a=0; a<nat; ++a )
    idxmap[a] = a;
  return idxmap;
}


class LList
{
public:

  LList();
  
  LList( int nat,
	 int const * pntarget,
	 std::vector<double> const & crds,
	 std::vector<int> const & idxmap );
  
  ~LList();

  int nat;
  int * ilist;
  int * numneigh;
  int * * firstneigh;
  deepmd::hpp::InputNlist * llist;
};


LList::LList()
  : nat(0),
    ilist(NULL),
    numneigh(NULL),
    firstneigh(NULL),
    llist(NULL)
{}


LList::~LList()
{
  if ( nat > 0 )
    {
      delete llist;
      delete [] ilist;
      delete [] numneigh;
      for ( int i=0; i<nat; ++i )
	delete [] firstneigh[i];
      delete [] firstneigh;
    };
}


LList::LList
( int inat,
  int const * pntarget,
  std::vector<double> const & crds,
  std::vector<int> const & idxmap )
  : nat(inat)
{
  double qmcut2 = ml::QMCUT*ml::QMCUT;
  std::vector< std::vector<int> > nlist( nat );
  for ( int i=1; i<*pntarget; ++i )
    for ( int j=0; j<i; ++j )
      {
	nlist[ idxmap[i] ].push_back( idxmap[j] );
	nlist[ idxmap[j] ].push_back( idxmap[i] );
      }
  for ( int i=0; i<*pntarget; ++i )
    for ( int j=*pntarget; j<nat; ++j )
      {
	double dx = crds[0+i*3]-crds[0+j*3];
	double dy = crds[1+i*3]-crds[1+j*3];
	double dz = crds[2+i*3]-crds[2+j*3];
	double r2 = dx*dx+dy*dy+dz*dz;
	if ( r2 < qmcut2 )
	  {
	    nlist[ idxmap[i] ].push_back( idxmap[j] );
	    nlist[ idxmap[j] ].push_back( idxmap[i] );
	  }
      }

  ilist = new int [nat];
  for ( int i=0; i<nat; ++i )
    ilist[i]=i;
  numneigh = new int [nat];
  for ( int i=0; i<nat; ++i )
    numneigh[i] = nlist[i].size();
  firstneigh = new int * [nat];
  for ( int i=0; i<nat; ++i )
    {
      firstneigh[i] = new int [numneigh[i]];
      for ( int j=0; j<numneigh[i]; ++j )
	firstneigh[i][j] = nlist[i][j];
    }
  llist = new deepmd::hpp::InputNlist( nat, ilist, numneigh, firstneigh );
}




static std::vector<std::string> get_typemap( std::string & orcstr )//std::string model )
{

  std::vector<std::string> names;
  std::string name;
  std::istringstream iss(orcstr.c_str());
  int aid=0;
  while( iss >> name )
    {
      names.push_back(name);
      //std::cout << aid << " " << names[aid] << "\n";
      aid++;
    };
  return names;
}




static double anint( double const x )
{
  return ((x)>0? std::floor((x)+0.5) : std::ceil((x)-0.5));
}


static void wrapcrds( double * R, double const * ucell, double const * recip )
{
  double frac[3] = { anint(R[0] * recip[0] + R[1] * recip[1] + R[2] * recip[2]),
		     anint(R[0] * recip[3] + R[1] * recip[4] + R[2] * recip[5]),
		     anint(R[0] * recip[6] + R[1] * recip[7] + R[2] * recip[8]) };
  R[0] -= frac[0] * ucell[0] + frac[1] * ucell[3] + frac[2] * ucell[6];
  R[1] -= frac[0] * ucell[1] + frac[1] * ucell[4] + frac[2] * ucell[7];
  R[2] -= frac[0] * ucell[2] + frac[1] * ucell[5] + frac[2] * ucell[8];
}


// static void wrapcrds( double * R, double const * ucell, double const * recip )
// {
//   double frac[3] = { anint(R[0] * recip[0] + R[1] * recip[3] + R[2] * recip[6]),
// 		     anint(R[0] * recip[1] + R[1] * recip[4] + R[2] * recip[7]),
// 		     anint(R[0] * recip[2] + R[1] * recip[5] + R[2] * recip[8]) };
//   R[0] -= frac[0] * ucell[0] + frac[1] * ucell[3] + frac[2] * ucell[6];
//   R[1] -= frac[0] * ucell[1] + frac[1] * ucell[4] + frac[2] * ucell[7];
//   R[2] -= frac[0] * ucell[2] + frac[1] * ucell[5] + frac[2] * ucell[8];
// }




extern "C"
{

  void new_dprc_( int const * ninter,
		char ** interfile,
		int const * nintra,
		char ** intrafile,
		int const * ntb,
		int const * calcavg //,
		//int const * usemmresnums
#ifdef WITH_MPI
		,int const * icommsander
#endif
		);

  void cpt_dprc_( double const * glb_crds,
		double * ml_ene,
		double * glb_frcs,
		double * maxstd,
		int const * pntarget,
		int const * pnat,
		int const * idxs,
		int const * residxs,
		char const ** atomnames,
		double const * ucell,
		int const * calcstd,
		int const * irank,
		int const * nrank );

  void img_dprc_( int const * pnat,
		double * crd,
		double const * ucell,
		double const * recip,
		int const * pnquant,
		int const * iqmatoms,
		double const * pqmcut,
		int * imm );
  
}








class glbmleval
{
public:
  
  glbmleval( std::vector<std::string> const & model_filename,
	     int const * ntb,
	     bool const calcavg,
	     bool interactor );
  
  
  void cpt_dprc( double const * glb_crds,
	       double * ml_ene,
	       double * glb_frcs,
	       double * maxstd,
	       int const * pntarget,
	       int const * pnat,
	       int const * idxs,
	       int const * residxs,
	       char const ** atomnames,
	       double const * ucell,
	       bool calcstd,
	       int const irank,
	       int const nrank );
  
private:
  
  bool interactor;
  int ntb;
  bool calcavg;
  int nmodel;
  std::vector<std::string> typemap;
  std::shared_ptr<deepmd::hpp::DeepPot> nnp_inter;
  std::shared_ptr<deepmd::hpp::DeepPotModelDevi> nnp_inter_devi;

  std::vector<double> virial;
  std::vector<double> box;
  
  std::vector<std::string> missing_atoms;
};








glbmleval::glbmleval
( std::vector< std::string > const & model_filename,
  int const * pntb,
  bool const icalcavg,
  bool inter )
  : interactor(inter),
    ntb( *pntb ),
    calcavg( icalcavg ),
    nmodel( model_filename.size() ),
    virial(9,0.),
    box(9,0.)
{

  if ( model_filename.size() == (std::size_t)1 )
    {
      //std::cout << "Model: " << model_filename[0] << std::endl;
      
      //nnp_inter.reset( new NNPInter(model_filename[0]) );
      nnp_inter.reset( new deepmd::hpp::DeepPot(model_filename[0]) );
      std::string orcstr;
      nnp_inter->get_type_map(orcstr);
      typemap = get_typemap( orcstr );
    }
  else if ( model_filename.size() > (std::size_t)1 )
    {
      //nnp_inter.reset( new NNPInter(model_filename[0]) );
      nnp_inter.reset( new deepmd::hpp::DeepPot(model_filename[0]) );
      std::string orcstr;
      nnp_inter->get_type_map(orcstr);
      typemap = get_typemap( orcstr );

      if ( calcavg )
	{
	  nnp_inter.reset();
	}
      
      //nnp_inter_devi.reset( new NNPInterModelDevi(model_filename) );
      nnp_inter_devi.reset( new deepmd::hpp::DeepPotModelDevi(model_filename) );
    }

  box[0] = BOXL;
  box[4] = BOXL;
  box[8] = BOXL;
  //nnp_inter->print_summary("ML summary: ");
}





//static int GLBCNT = 0;


void glbmleval::cpt_dprc
( double const * glb_crds,
  double * ml_ene,
  double * glb_frcs,
  double * maxstd,
  int const * pntarget,
  int const * pnat,
  int const * idxs,
  int const * residxs,
  char const ** atomnames,
  double const * nblist__ucell,
  bool calcstd,
  int const irank,
  int const nrank )
{
  //GLBCNT++;
  
  //double const KCAL_PER_EV = 23.06;
  double const KCAL_PER_EV = 1./0.04336410390059322;

  *ml_ene = 0.;
  *maxstd = 0.;

  if ( ! nnp_inter_devi )
    {
      calcstd = false;
    }
  else if ( calcavg )
    {
      calcstd = true;
    }
  
  std::vector<double> crds;
  std::vector<double> frcs;
  std::vector<int> atom_idxs;
  std::vector<int> atom_types;
  std::vector<std::string> atom_names;
  std::vector<int> residxs_vec;
  int nat=0;

  int natrange = *pntarget;
  if ( interactor )
    natrange = *pnat;


//  for ( int i=0; i<natrange; ++i )
//  {
//    std::cout << std::setw(6) << i << std::setw(6) << idxs[i] << "  " << atomnames[i] << std::endl;
//  }


  // std::cout << "calcavg " << calcavg << " " << calcstd
  // 	    << " " << (nnp_inter ? 1 : 0)
  // 	    << " " << (nnp_inter_devi ? 1 : 0)
  // 	    << std::endl;
  

  //std::printf("%i %i\n",ml::mpirank,ml::mpisize);
  
  for ( int i=0; i<natrange; ++i )
    {
      std::string name( atomnames[i] );
      std::vector<std::string>::iterator p = std::find( typemap.begin(), typemap.end(), name );
      if ( p != typemap.end() )
	{
	  nat++;
	  atom_types.push_back( std::distance( typemap.begin(), p ) );
	  atom_idxs.push_back( idxs[i] );
	  atom_names.push_back( name );
	  residxs_vec.push_back( residxs[i] );
	  
	  //f90stream::cout << "Atom " << std::setw(6) << i+1 << " named " << std::setw(4) << name << " has ML ID " << std::setw(3) << atom_types.back() << std::endl;
	  
	  for ( int k=0; k<3; ++k )
	    {
	      crds.push_back( glb_crds[k+idxs[i]*3] );
	      frcs.push_back(0.);
	    };
	}
      else if ( i < *pntarget )
	{
	  f90stream::cout << "Could not find target atom named " << name
			  << " within ML parameter file containing" << std::endl;
	  for ( std::vector<std::string>::iterator pp = typemap.begin(),
		  pend = typemap.end();
		pp != pend; ++pp )
	    {
	      f90stream::cout << " " << *pp;
	    }
	  f90stream::cout << std::endl;
	  std::exit(1);
	}
      else
	{
	  bool already_found_missing = false;
	  for ( std::size_t imiss=0; imiss < missing_atoms.size(); ++imiss )
	    {
	      if ( missing_atoms[imiss] == name )
		{
		  already_found_missing = true;
		  break;
		}
	    }
	  if ( ! already_found_missing )
	    {
	      missing_atoms.push_back(name);
	      f90stream::cout << "WARNING ML POTENTIAL MISSING PARAMETERS FOR ATOM: "
			      << name << std::endl;
	    }
	}
      
    }

  if ( ntb > 0 )
    {
      for ( int i=0; i<9; ++i )
	box[i] = nblist__ucell[i];
    }
  else
    {
      std::fill(box.data(),box.data()+9,0.);
      box[0] = BOXL;
      box[4] = BOXL;
      box[8] = BOXL;
    }

  ///
  // For some reason, setting the box to something big
  // will conserve energy in NVE, but setting the box
  // to Amber's true simulation cell (often quite small)
  // will not conserve energy
  //
  std::fill(box.data(),box.data()+9,0.);
  box[0] = BOXL;
  box[4] = BOXL;
  box[8] = BOXL;


  
  // {
  //   std::stringstream fname;
  //   fname << "frame_" << std::setfill('0') << std::setw(4) << GLBCNT << ".xyz";
  //   std::ofstream pdb;
  //   pdb.open( fname.str().c_str() );
  //   pdb << atom_names.size() << "\n";
  //   pdb << "\n";
  //   for ( std::size_t a=0, n=atom_names.size(); a<n; ++a )
  //     {
  // 	pdb << std::setw(4) << atom_names[a];
  // 	for ( std::size_t k=0; k<3; ++k )
  // 	  pdb << std::fixed << std::setw(12) << std::setprecision(7)
  // 	      << crds[k+a*3];
  // 	pdb << "\n";
  //     };
  //   pdb.close();
  // }


  {
    // NOW CPU SECTION

#ifdef WITH_MPI
    if ( nrank > 1 )
      {
	//std::cout << "nat " << nat << std::endl;

	std::vector<int> buckets(nrank,0);
	for ( int a=0; a<nat; ++a )
	  buckets[ a % nrank ]++;
	int ifirst = 0;
	for ( int a=0; a<irank; ++a )
	  ifirst += buckets[a];
	int ilast = ifirst + buckets[irank];
	
	//for ( int a=0; a<nrank; ++a )
	//std::cout << "bucket " << a << " = " << buckets[a] << std::endl;
	//std::cout << "rank " << irank << " " << ifirst << " " << ilast << std::endl;
	
	std::vector<double> scrds(3*nat,0.);
	std::vector<int> stypes(nat,0);
	std::vector<int> sresidxs(nat,0);
	std::vector<int> old2new(nat,0);
	{
	  int i=0;
	  for ( int a=ifirst; a<ilast; ++a, ++i )
	    {
	      old2new[a] = i;
	      stypes[i] = atom_types[a];
	      sresidxs[i] = residxs_vec[a];
	      for ( int k=0; k<3; ++k )
		scrds[k+i*3] = crds[k+a*3];
	    };
	  for ( int a=0; a<ifirst; ++a, ++i )
	    {
	      old2new[a] = i;
	      stypes[i] = atom_types[a];
	      sresidxs[i] = residxs_vec[a];
	      for ( int k=0; k<3; ++k )
		scrds[k+i*3] = crds[k+a*3];
	    };
	  for ( int a=ilast; a<nat; ++a, ++i )
	    {
	      old2new[a] = i;
	      stypes[i] = atom_types[a];
	      sresidxs[i] = residxs_vec[a];
	      for ( int k=0; k<3; ++k )
		scrds[k+i*3] = crds[k+a*3];
	    };
	  //std::cout << "i exits with " << i << std::endl;
	}
	int nghost = nat - buckets[irank];
	LList nlist( nat, pntarget, crds, old2new );
	int ago=0;
	
	*maxstd = 0.;
	
	if ( nnp_inter and not calcstd )
	  {
	    if (nnp_inter->dim_aparam() > 0)
	      {
		std::vector<double> sresidxsdouble(sresidxs.begin(), sresidxs.end());
		nnp_inter->compute( *ml_ene, frcs, virial, scrds, stypes, box, nghost,*(nlist.llist),ago,
				    std::vector<double>(), sresidxsdouble);
	      }
	    else
	      {
		nnp_inter->compute( *ml_ene, frcs, virial, scrds, stypes, box, nghost,*(nlist.llist),ago );
	      }
	    //nnp_inter->compute( *ml_ene, frcs, virial, scrds, stypes, sresidxs, box, nghost,*(nlist.llist),ago );
	    
	    for ( int a=0; a<nat; ++a )
	      {
		int i = old2new[a];
		int j = atom_idxs[a];
		for ( int k=0; k<3; ++k )
		  glb_frcs[k+j*3] += KCAL_PER_EV * frcs[k+i*3];
	      };
	    
	  }
	else if ( nnp_inter_devi )
	  {
	    
	    std::vector<double> all_ene( nmodel, 0. );
	    std::vector< std::vector<double> > all_frc(nmodel);
	    std::vector< std::vector<double> > all_virial(nmodel);
	    for ( int imodel=0; imodel<nmodel; ++imodel )
	      {
		all_frc[imodel].assign( frcs.size(), 0. );
		all_virial[imodel].assign( virial.size(), 0. );
	      };
	    // TODO: SWAP THESE LINES
	    if (nnp_inter_devi->dim_aparam() > 0)
	      {
		std::vector<double> sresidxsdouble(sresidxs.begin(), sresidxs.end());
		nnp_inter_devi->compute(all_ene,all_frc,all_virial,scrds,stypes,box,nghost,*(nlist.llist),ago,
					std::vector<double>(), sresidxsdouble);
	      }
	    else
	      {
		nnp_inter_devi->compute(all_ene,all_frc,all_virial,scrds,stypes,box,nghost,*(nlist.llist),ago);
	      }
	    //nnp_inter_devi->compute(all_ene,all_frc,all_virial,scrds,stypes,sresidxs,box,nghost,*(nlist.llist),ago);
	    
	    int ncrd = frcs.size();
	    std::vector<double> tmp(ncrd);
	    for ( int imodel=0; imodel<nmodel; ++imodel )
	      {
		for ( int a=0; a<nat; ++a )
		  {
		    int b=old2new[a];
		    for ( int k=0; k<3; ++k )
		      tmp[k+a*3] = all_frc[imodel][k+b*3];
		  }
		all_frc[imodel] = tmp;
		std::fill(tmp.begin(),tmp.end(),0.);
		MPI_Allreduce( all_frc[imodel].data(), tmp.data(), ncrd, MPI_DOUBLE, MPI_SUM, ml::comm );
		all_frc[imodel] = tmp;
		
		MPI_Allreduce( &(all_ene[imodel]), tmp.data(), 1, MPI_DOUBLE, MPI_SUM, ml::comm );
		all_ene[imodel] = tmp[0];
	      }
	    get_avg_and_std(calcavg,all_ene,all_frc,*ml_ene,frcs,*maxstd);
	    if ( irank > 0 )
	      {
		std::fill(frcs.begin(),frcs.end(),0.);
		*ml_ene = 0.;
	      }
	    
	    for ( int a=0; a<nat; ++a )
	      {
		int i = old2new[a];
		int j = atom_idxs[a];
		for ( int k=0; k<3; ++k )
		  glb_frcs[k+j*3] += KCAL_PER_EV * frcs[k+i*3];
	      };
	    
	  }
	else
	  {
	    std::cerr << "sander/ml.cpp failed to calc ML correction "
		      << "nnp_inter parallel CPU section nrank>1" << std::endl;
	    std::exit(1);
	  }
	
	*ml_ene *= KCAL_PER_EV;
	*maxstd *= KCAL_PER_EV;
	
	
      } // nrank > 1
    else
      {
#endif
	// cpu nrank == 1

	std::vector<int> old2new( MakeIdxMap(nat) );
	LList nlist( nat, pntarget, crds, old2new );
	int ago = 0;

	
	if ( nnp_inter and not calcstd )
	  {
	    if (nnp_inter->dim_aparam() > 0)
	      {
		std::vector<double> sresidxsdouble(residxs_vec.begin(), residxs_vec.end());
		nnp_inter->compute( *ml_ene, frcs, virial, crds, atom_types, box, 0, *(nlist.llist),ago, std::vector<double>(), sresidxsdouble);
	      }
	    else
	      {
		nnp_inter->compute( *ml_ene, frcs, virial, crds, atom_types, box, 0, *(nlist.llist),ago );
	      }
	    //nnp_inter->compute( *ml_ene, frcs, virial, crds, atom_types, residxs_vec, box );
	    *maxstd = 0.;
	  }
	else if ( nnp_inter_devi )
	  {
	    //std::vector<int> old2new( MakeIdxMap(nat) );
	    //LList nlist( nat, pntarget, crds, old2new );
	    //int ago = 0;
	    
	    std::vector<double> all_ene( nmodel, 0. );
	    std::vector< std::vector<double> > all_frc(nmodel);
	    std::vector< std::vector<double> > all_virial(nmodel);
	    for ( int imodel=0; imodel<nmodel; ++imodel )
	      {
		all_frc[imodel].assign( frcs.size(), 0. );
		all_virial[imodel].assign( virial.size(), 0. );
	      };
	    // TODO: SWAP THESE LINES
      if (nnp_inter_devi->dim_aparam() > 0) {
        std::vector<double> sresidxsdouble(residxs_vec.begin(), residxs_vec.end());
        nnp_inter_devi->compute(all_ene,all_frc,all_virial,crds,atom_types,box,0,*(nlist.llist),ago, std::vector<double>(), sresidxsdouble);
      } else {
	      nnp_inter_devi->compute(all_ene,all_frc,all_virial,crds,atom_types,box,0,*(nlist.llist),ago);
      }
	    //nnp_inter_devi->compute(all_ene,all_frc,all_virial,crds,atom_types,residxs_vec,box,0,*(nlist.llist),ago);
	    get_avg_and_std(calcavg,all_ene,all_frc,*ml_ene,frcs,*maxstd);
	  }
	else
	  {
	    std::cerr << "sander/ml.cpp failed to calc ML correction "
		      << "nnp_inter serial CPU section nrank==0" << std::endl;
	    std::exit(1);
	  }
	
	*ml_ene *= KCAL_PER_EV;
	*maxstd *= KCAL_PER_EV;
	
	for ( int i=0; i<nat; ++i )
	  {
	    int j = atom_idxs[i];
	    for ( int k=0; k<3; ++k )
	      glb_frcs[k+j*3] += KCAL_PER_EV * frcs[k+i*3];
	  };
	
#ifdef WITH_MPI
      }
#endif
  }
}



/////////////////////////////////////



std::shared_ptr< glbmleval > glbmlinter;
std::shared_ptr< glbmleval > glbmlintra;


void new_dprc_( int const * ninter,
	      char ** interfile,
	      int const * nintra,
	      char ** intrafile,
	      int const * ntb,
	      int const * pcalcavg //,
	      //int const * usemmresnums
#ifdef WITH_MPI
	      , int const * icommsander
#endif
	      )
{
  std::vector<std::string> inters;
  std::vector<std::string> intras;
  for ( int i=0; i<*ninter; ++i )
    inters.push_back( std::string( interfile[i] ) );
  for ( int i=0; i<*nintra; ++i )
    intras.push_back( std::string( intrafile[i] ) );

#ifdef WITH_MPI
  ml::comm = MPI_Comm_f2c( *icommsander );
  MPI_Comm_rank(ml::comm, &ml::mpirank);
  MPI_Comm_size(ml::comm, &ml::mpisize); 
#endif

  if ( *ninter > 0 )
    {
	{
	  glbmlinter.reset( new glbmleval( inters, ntb, *pcalcavg, true ) );
	};
    };
  if ( *nintra > 0 )
    {
	{
	  glbmlintra.reset( new glbmleval( intras, ntb, *pcalcavg, false ) );
	}
    };
}


void cpt_dprc_( double const * glb_crds,
	      double * ml_ene,
	      double * glb_frcs,
	      double * maxstd,
	      int const * pntarget,
	      int const * pnat,
	      int const * idxs,
	      int const * residxs,
	      char const ** atomnames,
	      double const * ucell,
	      int const * pcalcstd,
	      int const * irank,
	      int const * nrank )
{

  double maxstd_intra = 0.;
  double maxstd_inter = 0.;
  
  *ml_ene=0.;
  if ( glbmlinter )
    {
      double e=0.;
      glbmlinter->cpt_dprc( glb_crds, &e, glb_frcs, &maxstd_inter, pntarget, pnat,
			  idxs, residxs, atomnames, ucell, *pcalcstd, *irank, *nrank );
      *ml_ene += e;
    }
  
  if ( *irank == 0 )
    {
      if ( glbmlintra )
	{
	  double e=0.;
	  glbmlintra->cpt_dprc( glb_crds, &e, glb_frcs, &maxstd_intra, pntarget, pnat,
			      idxs, residxs, atomnames, ucell, *pcalcstd, 0, 1 );
	  *ml_ene += e;
	}
    };
      
  *maxstd = std::max( maxstd_inter, maxstd_intra );
}


void img_dprc_
( int const * pnat,
  double * crd,
  double const * ucell,
  double const * recip,
  int const * pnquant,
  int const * iqmatoms,
  double const * pqmcut,
  int * imm )
{
  //return;
  //double const A = 0.529177249;
  
  ml::QMCUT = *pqmcut;
  int const nat = *pnat;
  int const nquant = *pnquant;
  double const qmcut = *pqmcut;
  double const qmcut2 = qmcut*qmcut;
  std::fill( imm, imm + nat, 0 );
  


  //double ocrd[3] = { crd[0+iqmatoms[0]*3], crd[1+iqmatoms[0]*3],crd[2+iqmatoms[0]*3] };
  
  double avg[3] = { 0.,0.,0. };


  // for ( int i=0; i < nquant; ++i )
  //   f90stream::cout << iqmatoms[i] << "\n";
  
  std::vector<double> Rsums(nquant,0.);
  for ( int i=0; i < nquant; ++i )
    {
      int const a = iqmatoms[i];
      for ( int j=i+1; j<nquant; ++j )
  	{
  	  int const b = iqmatoms[j];
  	  double R[3] = { crd[0+a*3] - crd[0+b*3],
  			  crd[1+a*3] - crd[1+b*3],
  			  crd[2+a*3] - crd[2+b*3] };
  	  wrapcrds( R, ucell, recip );
  	  double Rmag=std::sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
	  Rsums[i] += Rmag;
	  Rsums[j] += Rmag;
  	};
    };

  double rmin = 1.e+10;
  int imin = 0;
  for ( int i=0; i < nquant; ++i )
    {
      if ( Rsums[i] < rmin )
	{
	  imin = i;
	  rmin = Rsums[i];
	};
    }
  for ( int k=0; k<3; ++k )
    avg[k] = crd[k+iqmatoms[imin]*3];

  wrapcrds( avg, ucell, recip );

  double box[6] = { 1.e+6, -1.e+6,  1.e+6, -1.e+6,  1.e+6, -1.e+6 };

  for ( int i=0; i < nquant; ++i )
    {
      int const a = iqmatoms[i];
      double R[3] = { crd[0+a*3] - avg[0],
		      crd[1+a*3] - avg[1],
		      crd[2+a*3] - avg[2] };
      wrapcrds( R, ucell, recip );
      box[0] = std::min( box[0], R[0] );
      box[1] = std::max( box[1], R[0] );
      box[2] = std::min( box[2], R[1] );
      box[3] = std::max( box[3], R[1] );
      box[4] = std::min( box[4], R[2] );
      box[5] = std::max( box[5], R[2] );
    };

  box[0] += avg[0]-qmcut;
  box[1] += avg[0]+qmcut;
  box[2] += avg[1]-qmcut;
  box[3] += avg[1]+qmcut;
  box[4] += avg[2]-qmcut;
  box[5] += avg[2]+qmcut;
  
  avg[0] = 0.5 * ( box[0]+box[1] );
  avg[1] = 0.5 * ( box[2]+box[3] );
  avg[2] = 0.5 * ( box[4]+box[5] );

  box[0] -= avg[0];
  box[1] -= avg[0];
  box[2] -= avg[1];
  box[3] -= avg[1];
  box[4] -= avg[2];
  box[5] -= avg[2];

  for ( int a = 0; a < nat; ++a )
    {
      crd[0+a*3] -= avg[0];
      crd[1+a*3] -= avg[1];
      crd[2+a*3] -= avg[2];
      wrapcrds( crd+a*3, ucell, recip );
    };
  
  int const * iqmatoms_end = iqmatoms + nquant;
  int const * b_end = iqmatoms + nquant;
  for ( int a=0; a < nat; ++a )
    {
      if ( crd[0+a*3] > box[0] and crd[0+a*3] < box[1] and
      	   crd[1+a*3] > box[2] and crd[1+a*3] < box[3] and
      	   crd[2+a*3] > box[4] and crd[2+a*3] < box[5] )
  	{
  	  int const * p = std::find( iqmatoms, iqmatoms_end, a );
  	  if ( p == iqmatoms_end ) // "a" is not a qm nor link atom
  	    {
  	      for ( int const * b=iqmatoms; b < b_end; ++b )
  		{
  		  double R[3] = 
  		    { crd[0+a*3]-crd[0+*b*3],
  		      crd[1+a*3]-crd[1+*b*3],
  		      crd[2+a*3]-crd[2+*b*3] };

  		  wrapcrds( R, ucell, recip );

  		  double const r2 = R[0]*R[0] + R[1]*R[1] + R[2]*R[2];

  		  if ( r2 < qmcut2  )
  		    {
		      imm[a] = 1;
  		      break;
  		    };
  		};
  	    };
  	};
      
    };

  {
    int i = iqmatoms[0];
    //double d[3] = { ocrd[0] - crd[0+i*3],
    //		    ocrd[1] - crd[1+i*3],
    //		    ocrd[2] - crd[2+i*3] };

    double d[3] = { 0.5*BOXL - crd[0+i*3],
    		    0.5*BOXL - crd[1+i*3],
    		    0.5*BOXL - crd[2+i*3] };
    
    for ( int a=0; a<nat; ++a )
      for ( int k=0; k<3; ++k )
	crd[k+a*3] += d[k];
  }

  
}




#endif
