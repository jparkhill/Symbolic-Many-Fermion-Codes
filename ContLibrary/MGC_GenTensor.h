#ifndef NOMGC

////////////////////////////////////////////////////////////////////////////////
// John Parkhill 2008
////////////////////////////////////////////////////////////////////////////////
#ifndef _GenTensor_C_
#define _GenTensor_C_

#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include "stdio.h"

#include "MGC_DenseLayered.h"
#include "MGC_filtering.h"
#include "MGC_Dense.h"
#include "MGC_Sets.h" 
#include "MGC_TensorSymmetry.h"
#include "MGC_CoordRep.h" 

using namespace std; 

// Numbers below apply to amplitude tensors, integrals, in ccman. 
class OrbitalTypes
{
public: 
	// Orbitals (static+dynamic)
	std::vector<int> oc,vr;
	// Active Orbitals 
	std::vector<int> activeo;
	std::vector<int> activev;
	// (dynamic)
	std::vector<int> dyno; 
	std::vector<int> dynv; 	
	// Frozen (above must be invoked)
	std::vector<int> frzno; 
	std::vector<int> frznv; 		
	// Orbital Spins
	std::vector<int> ospin;
	std::vector<int> vspin;	
	// Pair Map if it exists. 
	map<int,std::vector<int> > pair_map; 	
	// If there is sparsity this map contains the pair index of each spin-orbital
	map<int,int> PairForIndexO; 
	map<int,int> PairForIndexV; 
	
	OrbitalTypes()
	{} 
	OrbitalTypes(std::vector<int>& ao, std::vector<int>& av, 
				 std::vector<int>& dyo, std::vector<int>& dyv, 
				 std::vector<int>& fo, std::vector<int>& fv,
				 std::vector<int>& os, std::vector<int>& vs, 
				 map<int,std::vector<int> >& pm): activeo(ao), activev(av), dyno(dyo), dynv(dyv),
	frzno(fo), frznv(fv), ospin(os), vspin(vs), pair_map(pm)
	{
		MakePairForIndex();
		return;
	}
	OrbitalTypes(std::vector<int>& o, std::vector<int>& v, 
				 std::vector<int>& ao, std::vector<int>& av, 
				 std::vector<int>& dyo, std::vector<int>& dyv, 
				 std::vector<int>& fo, std::vector<int>& fv,
				 std::vector<int>& os, std::vector<int>& vs, 
				 map<int,std::vector<int> >& pm): oc(o), vr(v), activeo(ao), activev(av), dyno(dyo), dynv(dyv),
	frzno(fo), frznv(fv), ospin(os), vspin(vs), pair_map(pm)
	{
		MakePairForIndex();
		return;
	}
	OrbitalTypes& operator=(const OrbitalTypes& other)
	{
		oc = other.oc; vr = other.vr; 
		activeo = other.activeo; activev = other.activev; 
		dyno = other.dyno; dynv = other.dynv; 
		frzno = other.frzno; frznv = other.frznv; 
		ospin=other.ospin; vspin = other.vspin; 
		pair_map = other.pair_map; 
		PairForIndexO = other.PairForIndexO; 
		PairForIndexV = other.PairForIndexV;
		return *this;
	}
	void clear()
	{
		oc.clear(); vr.clear();
		activeo.clear(); activev.clear();
		dyno.clear(); dynv.clear();
		frzno.clear(); frznv.clear();
		ospin.clear(); vspin.clear();
		pair_map.clear();
		PairForIndexO.clear();
		PairForIndexV.clear();
	}
	void MakePairForIndex()
	{
		PairForIndexO.clear(); 
		PairForIndexV.clear(); 
		for (map<int,std::vector<int> >::const_iterator iter = this->pair_map.begin(); iter!= this->pair_map.end(); iter++)
		{	
			for(int i = 0; i < 2; i++)
			{
				PairForIndexO[((iter->second)[i])] = iter->first; 
			}
			for(int i = 2; i < 4; i++)
			{
				PairForIndexV[((iter->second)[i])] = iter->first; 
			}
		}
		// if there are external do that mapping too. 
		if (dyno.size() || dynv.size())
		{
			std::vector<int>::const_iterator eiter; 
			for (eiter = dyno.begin(); eiter != dyno.end(); ++eiter)
			{
				PairForIndexO[*eiter] = -1;
			}
			for (eiter = dynv.begin(); eiter != dynv.end(); ++eiter)
			{
				PairForIndexV[*eiter] = -1;
			}			
		}
		return;
	}
	inline bool IfRestrOcc(const int& arg)
	{	
		return (std::find(activeo.begin(), activeo.end(), arg) == activeo.end());
	}
	inline bool IfRestrVirt(const int& arg)
	{
		return (std::find(activev.begin(), activev.end(), arg) == activev.end());
	}
	inline int NMO() const 
	{
		return (NAlpha()+NBeta()); 
	}
	inline int NBeta() const 
	{
		return (NBetaOcc()+NBetaVirt());
	}
	inline int NAlpha() const 
	{
		return (NAlphaOcc()+NAlphaVirt());
	}
	inline int NAlphaVirt() const
	{
		int tore=0; 
		for (std::vector<int>::const_iterator tmp = vspin.begin(); tmp != vspin.end(); ++tmp)
			if (*tmp == 0)
				++tore; 
		return tore; 
	}
	inline int NBetaVirt() const
	{
		int tore=0; 
		for (std::vector<int>::const_iterator tmp = vspin.begin(); tmp != vspin.end(); ++tmp)
			if (*tmp == 1)
				++tore; 
		return tore; 
	}	
	inline int NAlphaOcc() const
	{
		int tore=0; 
		for (std::vector<int>::const_iterator tmp = ospin.begin(); tmp != ospin.end(); ++tmp)
			if (*tmp == 0)
				++tore; 
		return tore; 
	}
	inline int NBetaOcc() const
	{
		int tore=0; 
		for (std::vector<int>::const_iterator tmp = ospin.begin(); tmp != ospin.end(); ++tmp)
			if (*tmp == 1)
				++tore; 
		return tore; 
	}	
	// args are c-matrix number of an orbital
	bool SameSpaceCMat(const int& p, const int& q, bool IfBeta=false)
	{
		return ( TypeOfCMat(p,IfBeta) == TypeOfCMat(q,IfBeta) );
	}
	// arg is c-matrix number of an orbital
	// result is an integer representing the type
	// 1 = inactive o, 2 = active o , 3 = active v , 4 = inactive virtual. 
	// CCman's Orbital numbering is as follows:  
	//  Occ:  0 |Rst. Alph.|Act. Alph.|Rst. Bet.|Act. Bet.| NOcc  -1
	// Virt:  0 |Act. Alph.|Rst. Alph.|Act. Bet.|Rst. Bet.| NVirt -1
	// Cmatrix Numbering is for (A and B separately)
	// (0) | R Occ | Active O | Active V | R Virt | (NMO-1)
	int TypeOfCMat(const int& arg, bool IfBeta = false ) const
	{
		int pt = arg; 
		if (!IfBeta)
		{
			if ( pt < NAlphaOcc() )
			{
				if ( std::find(activeo.begin(), activeo.end(), pt) != activeo.end() )
					return 2;
				else 
					return 1; 
			}
			else
			{
				pt -= NAlphaOcc();
				if ( std::find(activev.begin(), activev.end(), pt) != activev.end() )
					return 3;
				else 
					return 4; 				
			}
		}
		else
		{
			if ( pt < NBetaOcc() )
			{
				pt += NAlphaOcc(); 
				if ( std::find(activeo.begin(), activeo.end(), pt) != activeo.end() )
					return 2;
				else 
					return 1; 	
			}
			else 
			{
				pt -= NBetaOcc(); 
				pt += NAlphaVirt(); 
				if ( std::find(activev.begin(), activev.end(), pt) != activev.end() )
					return 3;
				else 
					return 4; 								
			}
		}
	}
	// Converts an index into it's absolute number: ie: 
	// 0, AlphaOcc1,...,AlphaVirt1,...N ... BetaOcc1 ...   
	std::vector<int> Untype(const std::vector<int>& arg,const std::vector<int>& otype) const
	{
		std::vector<int> tore(arg.size()); 
		if (arg.size() > otype.size())
		{
			cout << "Size mismatch in Untype: " << arg.size() << " " << otype.size() << endl;  
			return tore; 
		}
		for (std::vector<int>::const_iterator iter = arg.begin(); iter != arg.end(); ++iter)
			tore[iter-arg.begin()] = Untype(*iter,otype[iter-arg.begin()]);
		return tore; 
	}
	void Untype(const std::vector<int>& arg,const std::vector<int>& otype, std::vector<int>& result) const
	{
		if (arg.size() > otype.size())
		{
			cout << "Size mismatch in Untype: " << arg.size() << " " << otype.size() << endl;  
		}
		for (std::vector<int>::const_iterator iter = arg.begin(); iter != arg.end(); ++iter)
			result[iter-arg.begin()] = Untype(*iter,otype[iter-arg.begin()]);
		return ; 
	}
	bool UntypedO(const int& arg)
	{
		int tmp(arg); 
		if (arg < NAlpha())
		{
			if (arg < NAlphaOcc())
				return true; 
			else 
				return false; 
		}
		else 
		{
			if (arg < (NAlpha()+NBetaOcc()))
				return true; 
			else 
				return false; 			
		}
	}
	int Untype(const int& arg, const bool OorV) const
	{
		int tore;
		if (OorV)
		{
			// virt
			if (vspin[arg] == 0)
				tore = arg+NAlphaOcc();// alpha
			else 
				tore = arg+(NAlphaOcc()+NBetaOcc()); // beta NAlphaVirt Already in there. 
		}
		else 
		{
			if (ospin[arg] == 0)
				tore = arg;// alpha
			else 
				tore = arg+NAlphaVirt(); // beta (NAlphaOcc already in there)
		}
		return tore; 
	}
	int Retype(const int& arg)
	{
		int tore(arg); 
		if (tore < NAlpha())
		{
			if (tore < NAlphaOcc())
				return tore;
			else 
				tore -= NAlphaOcc(); 		
		}
		else 
		{
			tore -= NAlpha(); 
			if (tore < NBetaOcc())
				return (tore+NAlphaOcc());
			else 
				tore -= (NBetaOcc() - NAlphaVirt());
		}
		return tore; 
	}
	void Retype(std::vector<int>& arg,const int NtoRetype)
	{
		int NTR=0;
		for (std::vector<int>::iterator iter = arg.begin(); iter != arg.end(); ++iter)
		{
			if (NTR >= NtoRetype)
				break; 
			*iter = Retype(*iter);
			++NTR;
		}
	}	
	void Retype(std::vector<int>& arg)
	{
		std::vector<int> tore(arg); 
		for (std::vector<int>::iterator iter = tore.begin(); iter != tore.end(); ++iter)
			*iter = Retype(arg[iter-tore.begin()]);
		std::copy(tore.begin(),tore.end(),arg.begin()); 
	}
	std::vector<int> Retype(const std::vector<int>& arg)
	{
		std::vector<int> tore(arg); 
		for (std::vector<int>::iterator iter = tore.begin(); iter != tore.end(); ++iter)
			*iter = Retype(arg[iter-tore.begin()]);
		return tore; 
	}
	
};

/*
 Defines what functions are required to use Contraction.C 
 */

class GenTensor
{
public:
	int dimensions;
	
	//Should be made the way this information is kept everywhere. 
	OrbitalTypes MyTypes; 
	
	TensorSymmetry MySym;
	Filters MyFilter; 
	// Variables that control filtering & element addition. 
	bool if_paired; //Pairing enforced  
	int num_pairs; //Number of local pairs
	int if_dyn; // Amplitude filter for (+SD) model. 1 = Only into external 2 = Into any space. 
	bool if_frzn; // There are frozen orbitals.  
	
	// indicates what this tensor represents...
	bool IsInter; //Is this Tensor an intermediate? (If so del() does more)
	bool IsIntegral;
	bool IsAmp; 
	
	std::vector<int> o_type; // 0 = occ, 1 = virt, -1 = any. 
	std::vector<int> BraketType; // 1 Up // -1 Down 
	
	GenTensor* MyClone;// Useful for debugging. 	
	
	GenTensor() : if_paired(false), num_pairs(0), if_dyn(false), if_frzn(false), MyFilter(0), 
	IsInter(false), IsAmp(false), IsIntegral(false), dimensions(0)
	{
	}
	
	GenTensor(const int rank) : if_paired(false), num_pairs(0), dimensions(rank) , o_type(rank,0), MySym(rank),
	if_dyn(false), if_frzn(false), MyFilter(rank), IsInter(false), IsAmp(false), IsIntegral(false), BraketType(rank,1)
	{
	}
	
	GenTensor(const GenTensor& a_CopyOf) : o_type(a_CopyOf.o_type), dimensions(a_CopyOf.dimensions), MySym(a_CopyOf.MySym), MyFilter(a_CopyOf.MyFilter), 
	if_paired(a_CopyOf.if_paired), num_pairs(a_CopyOf.num_pairs), if_dyn(false), if_frzn(false), IsInter(a_CopyOf.IsInter), 
	IsAmp(a_CopyOf.IsAmp), IsIntegral(a_CopyOf.IsIntegral), BraketType(a_CopyOf.BraketType), MyTypes(a_CopyOf.MyTypes)
	{
	}
	
	virtual ~GenTensor()
	{
	}
	
	void TakeOrbitalInformation(const GenTensor& Source)
	{
		if_paired = Source.if_paired; 
		num_pairs = Source.num_pairs; 
		if_dyn = Source.if_dyn; 
		if_frzn = Source.if_frzn; 
		IsInter = Source.IsInter; 
		IsIntegral = Source.IsIntegral; 
		IsAmp = Source.IsAmp; 
		o_type = Source.o_type;
		BraketType = Source.BraketType; 
		MyTypes = Source.MyTypes; 
	}
	
	virtual GenTensor* clone() const = 0; 
	// Useful For Debuggin. 
	virtual std::vector<int> FirstInequality(GenTensor&) = 0 ;
	virtual bool operator==(GenTensor& other) = 0; 
	virtual const int Rank() const = 0; 
	virtual double l2norm() const = 0; 
	virtual int size() const = 0; 
	inline virtual bool obeys_filters(const std::vector<int>&) = 0 ; 
	virtual void ReStruc(double thresha = 0.0) = 0 ;
	virtual void print(int code = 10) const = 0 ; 
	virtual void incr(const std::vector<int>& a_size, const double& a_value) = 0; 
	virtual double val(const std::vector<int>& ) = 0 ; 
	virtual void set2zero() = 0 ; 
	virtual void ComplementaryDims(std::vector<int>& arg) = 0; 
	
	//Functions that provide simple raster ordering. 
	virtual bool nonzero(const std::vector<int>& ) const = 0; 
	virtual bool begin() = 0; 
	virtual double val() = 0 ; 
	virtual void set(const double& arg) = 0 ; 
	virtual std::vector<int> ind() = 0 ;
	virtual bool I() = 0 ;
	
	virtual DenseCoordRep* GetDCR() = 0 ; 
	virtual DenseLayeredRep* GetDLR(const std::vector<char*>& BlockingList = std::vector<char*>()) = 0; 
	virtual DenseLayeredRep* GetDLR(const OrbitalTypes& pio, const std::vector<char*>& BlockingList = std::vector<char*>()) = 0; 
};

#endif 
#endif 
