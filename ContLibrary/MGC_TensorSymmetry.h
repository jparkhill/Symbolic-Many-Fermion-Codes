#ifndef NOMGC

#ifndef tsm_C
#define tsm_C
/*
	This Class can only handle sets of indices which are anti-symmetric to each other. 
	In particular (NOT-ANTI) symmetries which sometimes exist on things like VOVOV 
	are ignored in my code in general (c.f.) BraKetHermitize() in MBTensor.C 
	
	Symmetry is Assigned with this process: 
	Sets of symmetrical dimensions => Groups on those dimensions => A complete group. 
	
	If something has O/V amplitude symmetry this is done with separate PGs for O,V & 
	sorting is made more efficient.  
	If it has random intermediate symmetry then sorting needs the dset and is done solely 
	through Ogroup. 
	
	JAP 2009. 
*/

#include <climits>
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include "MGC_StoreRecall.h"
#include "MGC_Sets.h"
#include "MGC_Permutation.h"

using namespace std; 

//Container for a tensor's permutational index symmetries. 
class TensorSymmetry
{
	public : 
	bool empty; 
	bool OneGroup; //Uses ppg. 
	int rank; //Rank of representing tensor
	std::map<int, std::vector<int> > dsym; // For each dimension, those dimensions which are (anti)symmetric to this one. 
	std::vector<std::vector<int> > dset; // (exclusive) Sets of dimensions symmetric to one another. 
	//These can be determined from above. 
	PermutationGroup OccPermGroup,VirtPermGroup; 
	ProductPermutationGroup ppg; 
	
	//Complete Symmetry for first and half subset of dimensions
	inline TensorSymmetry(int rank_arg = 0) : rank(rank_arg) , empty(true), OneGroup(true)
	{
		return; 
	}

	inline TensorSymmetry(const TensorSymmetry& other) : OccPermGroup(other.OccPermGroup) , VirtPermGroup(other.VirtPermGroup) , dsym(other.dsym) , dset(other.dset) , 
		rank(other.rank), empty(other.empty), OneGroup(other.OneGroup), ppg(other.ppg)
	{
		return; 
	}

	inline void clear()
	{
		empty = true; 
		OneGroup = true; 
		dsym.clear(); 
		dset.clear(); 
		OccPermGroup.clear(); 
		VirtPermGroup.clear(); 
		ppg.clear(); 
	} 
	/*
	std::vector<NotJustAny> MakeSummary() const
	{
		std::vector<boost::any> tore; 
		tore.push_back(empty);
		tore.push_back(OneGroup);		
		tore.push_back(rank);		
		tore.push_back(dsym);
		tore.push_back(dset);
	}
	void AssignFromSummary(std::vector<NotJustAny>& Arg) 
	{
		clear(); 
		if (Arg.size() < 5)
			cout << "Gibberish Summary";
		empty = boost::any_cast<bool>(Arg[0]);
		OneGroup = boost::any_cast<bool>(Arg[1]);
		rank = boost::any_cast<int>(Arg[2]);
		dsym = boost::any_cast<std::map<int, std::vector<int> > >(Arg[3]);
		dset = boost::any_cast<	std::vector<std::vector<int> > >(Arg[4]);
		if (OneGroup)
			AssignPG(dsym);
		else 
			AssignCompletePG();
	}
	*/
	inline void CopyFrom(const TensorSymmetry& other)
	{
		clear(); 
		if (!other.empty)
		{
			OccPermGroup.CopyFrom(other.OccPermGroup);
			VirtPermGroup.CopyFrom(other.VirtPermGroup);		
			if (other.ppg.rank)
				ppg.CopyFrom(other.ppg);
		}
		std::copy(other.dsym.begin(), other.dsym.end() , std::inserter(dsym,dsym.begin())); 
		dset.assign(other.dset.begin(), other.dset.end()); 
		rank = other.rank; 
		empty = other.empty;
		OneGroup = other.OneGroup; 
		return; 
	}
	bool operator==(const TensorSymmetry& other) const
	{
		if (empty != other.empty)
		{
			cout << "EmptyDiff" << endl; 
			return false; 
		}
		else if (dsym != other.dsym)
		{
			cout << "!=Dsym" << endl; 
			return false; 
		}
		else if (dset != other.dset)
		{
			cout << "!=Dsym" << endl; 
			return false; 
		}		
		else if (rank != other.rank)
		{
			cout << "!=rank" << endl; 
			return false; 
		}
		else if (OneGroup != other.OneGroup)
		{
			cout << "!=OneGroup" << endl; 
			return false; 
		}
		else if (!(OccPermGroup == other.OccPermGroup) || !(VirtPermGroup == other.VirtPermGroup) || !(ppg == other.ppg) )
		{
			cout << "!=PGroup" << endl; 
			return false; 
		}
		return true; 
	}
	
	inline long int order() const 
	{
		if (!empty && !OneGroup)
			return (OccPermGroup.size()*VirtPermGroup.size()); 
		else if (!empty && OneGroup)
			return (ppg.size()); 
		else 
			return 1; 
	} 
	
	inline void AssignPG(const std::map<int,std::vector<int> >& arg_dsym)
	{
		dsym.clear(); 
		dset.clear(); 
		std::copy(arg_dsym.begin(), arg_dsym.end() , std::inserter(dsym,dsym.begin())); 
		empty = true; 
		for (int i = 0 ; i < rank ; ++i)
		{
			if (dsym[i].size() != 0)
				empty = false; 
		}
		if (!empty)
		{
			//rank = dsym.size();
			dset = DsymToSets(dsym);
			ppg.clear(); 
			ppg = ProductPermutationGroup(rank, dset);
		}
		else 
		{
			ppg.clear(); 	
			empty = true;
		}
		OneGroup = true; 
		OccPermGroup.clear(); 
		VirtPermGroup.clear(); 
		return; 
	}
	
	void AssignCompletePG()
	{
		dsym.clear(); 
		dset.clear(); 
		for (int dim1 = 0 ; dim1 < rank ; dim1++)
		{
			for (int dim2 = dim1 + 1; dim2 < rank; dim2++)
			{
				if ((dim1 >= rank/2) || (dim2 >= rank/2))
					continue;
				(dsym[dim1]).push_back(dim2);
				(dsym[dim2]).push_back(dim1);
			}
		}
		for (int dim1 = 0 ; dim1 < rank ; dim1++)
		{
			for (int dim2 = dim1 + 1; dim2 < rank; dim2++)
			{
				if ((dim1 < rank/2) || (dim2 < rank/2))
					continue;
				(dsym[dim1]).push_back(dim2);
				(dsym[dim2]).push_back(dim1);
			}
		}
		this->dset = DsymToSets(this->dsym);
		OccPermGroup.GiveRank(rank/2);
		VirtPermGroup.GiveRank(rank/2);
		ppg.clear(); 
		OneGroup = false; 
		empty = false; 
	}

	// Special purpose sorting routine.
	// Pass in a list of contracted dimensions/
	// this gets them sorted and generates the corresponding permutation and sign. 
	// (See routine below)
	void SortCdim(std::vector<int>& ToSort)
	{ 
		std::vector<int> Tmp(rank,-1);
		std::vector<int>::iterator viter; 
		if (OneGroup)
		{
			for (viter = ToSort.begin(); viter != ToSort.end(); ++viter)
				Tmp[*viter] = (viter - ToSort.begin());
			int HoldToSwap = 0 ; 
			bool swapped = true; 
			std::vector<int>::iterator i,j;
			int iint = 0;
			int jint = 0;
			std::vector<int>* ISym; 
			while (swapped)
			{
				swapped = false; 
				for (i = Tmp.begin(); i != Tmp.end(); ++i)
				{
					iint = (i - Tmp.begin());
					// Shorthand
					for (j = i+1; j != Tmp.end(); ++j)
					{
						jint = (j - Tmp.begin());
						if (*i < 0) 
							break; 
						if (*j < 0)
							continue; 
						// does dsym allow this swap? 
						if (*i > *j)
						{
							if (std::find(dsym[iint].begin(), dsym[iint].end(), jint) != dsym[iint].end() )
							{
								swapped = true;
								std::iter_swap(i,j); 
							}
						}
						
					}					
				}
			}
		} 
		else
		{
			// Make sure bra precedes ket. 
			std::vector<int>::iterator i,j;
			bool swapped = true; 
			while (swapped)
			{
				swapped = false; 
				for (i=ToSort.begin(); i != ToSort.end(); ++i)
				{
					for (j=i+1 ; j != ToSort.end() ; ++j)
					{
						if (*i > rank/2 && *j < rank/2)
						{
							std::iter_swap(i,j);
							swapped = true; 
						}
					}				
				}
			}			
			for (viter = ToSort.begin(); viter != ToSort.end(); ++viter)
				Tmp[*viter] = (viter - ToSort.begin());			
			std::vector<int> inputo(rank/2), inputv(rank/2);
			std::vector<int> toreo(rank/2), torev(rank/2);
			for (int tmp = 0 ; tmp < rank/2; ++tmp)
			{
				toreo[tmp] = tmp;
				torev[tmp] = tmp; 
				inputo[tmp] = Tmp[tmp];
				inputv[tmp] = Tmp[tmp+(rank/2)];
			}
			int HoldToSwap = 0 ; 
			swapped = true; 
			while (swapped)
			{
				swapped = false; 
				for (i = inputo.begin(); i != inputo.end(); ++i)
				{
					// Shorthand
					for (j = i+1; j != inputo.end(); ++j)
					{
						if (*i < 0) 
							break; 
						if (*j < 0)
							continue; 
						// does dsym allow this swap? 
						if (*i > *j)
						{
							swapped = true;
							std::iter_swap(i,j); 
						}						
					}					
				}
			}
			swapped = true; 
			while (swapped)
			{
				swapped = false; 
				for (i = inputv.begin(); i != inputv.end(); ++i)
				{
					// Shorthand
					for (j = i+1; j != inputv.end(); ++j)
					{
						if (*i < 0) 
							break; 
						if (*j < 0)
							continue; 
						// does dsym allow this swap? 
						if (*i > *j)
						{
							swapped = true;
							std::iter_swap(i,j); 
						}						
					}					
				}
			}
			Tmp = Join(inputo,inputv); 
		}
		ToSort = OrderOfDimensionsToDimensionsInOrder(Tmp);
		return; 
	}	
	
	// Permutes CPair every which way, that sorts the contracted dimensions. 
	// and incorporates the appropriate factor. 
	void ExCdPreservingOrder(const std::vector<int>& MyCDims, 
							 std::vector<std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > >& MyCmap, 
							 int print_lvl = 0)
	{
		int CdimFactor = 1; 		
		std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > to_push;		
		std::vector<int> MyContDims(rank, -1);
		std::vector<int>::const_iterator viter,viter1; 		
		for (viter1 = MyCDims.begin(); viter1 != MyCDims.end(); ++viter1)
			MyContDims[*viter1] = (viter1-MyCDims.begin());
		if (print_lvl)
		{
			cout << " MyCDims: " ;
			printvec(MyCDims,true);
			cout << " MyContDims: " ;
			printvec(MyContDims,true);
		}
		
		CdimFactor = ObtainCdimFactor(MyContDims, print_lvl);
		if (print_lvl)
			cout << "ExCdPreservingOrder -- Factor: " << CdimFactor << endl; 
		//Construct 0-rank range which will be reused throughout. 
		std::vector<int> tmp(rank), rang(rank), tmp2(rank), TstVec(rank);
		for (int dim = 0 ; dim < rank ; dim++)
			rang[dim] = dim;
		std::vector<int> T3(MyCDims.size());
		std::set< std::vector<int> > TmpSet;  // To reduce stupidity ;P 
		bool SortedInCdims = true;
		if (empty || rank == 2)
		{
			if (print_lvl)
				cout << " ---- EmptyOrOneGroup in ExCdPreservingOrder  ---- " << endl; 
			to_push = std::make_pair(MyCDims,std::make_pair(rang,rang));
			MyCmap.push_back(to_push);
			MyCmap.back().second.first.push_back(1); 
			MyCmap.back().second.second.push_back(1); 
			return; 
		}
		else if (OneGroup)
		{
			ppg.Begin();
			while(true)
			{
				tmp = rang; 
				MutVectorShiftInPlace(tmp, ppg.Current() ,0);	
				tmp2 = SuckDims(MyCDims,tmp);
				if (MyCDims.size() > 1)
				{
					TstVec.assign(rank, -1); 
					for (viter = tmp2.begin(); viter != tmp2.end(); ++viter)
						TstVec[*viter] = (viter - tmp2.begin());
					SortedInCdims = ObeysMyOrdering(TstVec); 
				}
				else 
					SortedInCdims = true;	
				if (SortedInCdims)
				{
					if (TmpSet.count(tmp2) == 0)
					{
						to_push = std::make_pair(tmp2, std::make_pair(ppg.Current(),rang));
						to_push.second.second.push_back(1); // Second Permutation is just the identity. 
						to_push.second.first.back() *= CdimFactor;
						MyCmap.push_back(to_push);						
						if (print_lvl)
						{
							cout << "ExCdPres Order - 1grp Cmap I'm Pushing: " << endl; 
							printvec(to_push.first);
							printvec(to_push.second.first);
							printvec(to_push.second.second);
						}														
						TmpSet.insert(tmp2); 
					}
				}
				if (!ppg.Next())
					break; 
			}
			return; 
		}
		else 
		{
			int HlfRnk=rank/2;
			OccPermGroup.Begin();
			while(true)
			{
				VirtPermGroup.Begin();
				while(true)
				{
					for (int d = 0;d<MyCDims.size();++d)
					{
						if (MyCDims[d] < HlfRnk)
							T3[d]=OccPermGroup.Current(MyCDims[d]);
						else
							T3[d]=(VirtPermGroup.Current(MyCDims[d]-HlfRnk))+HlfRnk;
					}					
					if (MyCDims.size() > 1)
					{
						TstVec.assign(rank, -1); 
						for (viter = T3.begin(); viter != T3.end(); ++viter)
							TstVec[*viter] = (viter - T3.begin());
						SortedInCdims = ObeysMyOrdering(TstVec); 
					}
					else 
						SortedInCdims = true; 
					if (SortedInCdims)
					{					
						//Figure out if this presents a new Contraction. 
						if (TmpSet.count(T3) == 0)
						{
							to_push = std::make_pair(T3, std::make_pair(OccPermGroup.Current(),VirtPermGroup.Current()));
							MyCmap.push_back(to_push);
							MyCmap.back().second.first.back() *= CdimFactor;
							if (print_lvl)
							{
								cout << "ExCdPres Order - Cmap I'm Pushing: " << endl; 
								printvec(to_push.first);
								printvec(to_push.second.first);
								printvec(to_push.second.second);
							}														
							TmpSet.insert(T3); 
						}
					}
					if (!VirtPermGroup.Next())
						break; 
				}				
				if (!OccPermGroup.Next())
					break; 
			}
			return; 
		}
	}	
	
	// On the tensor of a contraction whose symmetries
	// contain the other's run this routine. 
	// it expands cpair for the symmetries on the contracted dimensions of this 
	// tensor which have not yet been accounted for. 
	void ExpandCpair( const std::vector<int>& MyCdims_arg, const std::vector<int>& OtherCdims, TensorSymmetry& Other,
					 std::vector<std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > >& MyCmap, 
					 std::vector<std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > >& OtherCmap, int print_lvl = 0)
	{
		MyCmap.clear(); 
		OtherCmap.clear(); 
		Other.ExCdPreservingOrder(OtherCdims, OtherCmap, print_lvl); 
		std::vector<int> SortedOtherCdims(OtherCdims);
		Other.SortCdim(SortedOtherCdims);
		if (print_lvl)
		{
			cout << "SortedOtherCdims " << endl; 
			printvec(SortedOtherCdims); 
		}			
		// OtherCmap.push_back(OtherCPair); 
		if (OneGroup)
		{
			cout << "ExpandCpair Being Called on a OneGroup Tensor" << endl;  
			system("sleep 100");
		}
		// Make sure occupieds precede virtuals. 
		std::vector<int> MyCdims(MyCdims_arg);
		std::vector<int>::iterator i,j;
		bool swapped = true; 
		while (swapped)
		{
			for (i=MyCdims.begin(); i != MyCdims.end(); ++i)
			{
				swapped = false; 
				for (j=i+1 ; j != MyCdims.end() ; ++j)
				{
					if (*i > rank/2 && *j < rank/2)
					{
						std::iter_swap(i,j);
						swapped = true; 
					}
				}				
			}
		}
		// Construct a vector to see if this obeys other tensor's symmetry. 
		std::vector<int> STest(Other.rank, -1);
		// Construct shifted permutation groups		
		PermutationGroup TheOGroup(rank/2);
		PermutationGroup TheVGroup(rank/2);
		// Construct the Cdims. 
		std::vector<int> T3(MyCdims.size()); 
		int HlfRnk = rank/2;
		bool ObeysOtherSymmetries = true; 
		std::set< std::vector<int> > TmpSet; // No repeat contractions. 
		std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > to_push;		
		TheOGroup.Begin();
		while(true)
		{
			TheVGroup.Begin();
			while(true)
			{
				for (int d = 0;d<MyCdims.size();++d)
				{
					if (MyCdims[d] < HlfRnk)
						T3[d]=TheOGroup.Current(MyCdims[d]);
					else
						T3[d]=(TheVGroup.Current(MyCdims[d]-HlfRnk))+HlfRnk;
				}
				// are in the same order... ?
				ObeysOtherSymmetries = true; 
				if (T3.size() > 1)
				{
					STest.assign(Other.rank, -1); 
					for (int i = 0; i < T3.size() ; ++i)
						STest[SortedOtherCdims[i]] = T3[i];  
					ObeysOtherSymmetries = Other.ObeysMyOrdering(STest); 
				}
				else 
					ObeysOtherSymmetries = true; 
				if (ObeysOtherSymmetries)
				{
					if (TmpSet.count(T3) == 0)
					{
						to_push = std::make_pair(T3, std::make_pair(TheOGroup.Current(), TheVGroup.Current()));
						if (print_lvl)
						{
							cout << "Cmap I'm Pushing: " << endl; 
							printvec(T3);
							printvec(TheOGroup.Current());
							printvec(TheVGroup.Current());
						}
						MyCmap.push_back(to_push);
						TmpSet.insert(T3); 
					}
				}				
				if (!TheVGroup.Next())
					break; 
			}				
			if (!TheOGroup.Next())
				break; 
		}		
		return; 
	}	
	
	// Make the Perm-Exploded equivalent of Cmap. 
	// That-is all contractions which are equivalent to this by symmetry
	// Yet no longer equivalent since we represent only sorted indices. 
	std::vector<std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > > ExplodeCdims(std::vector<int>& dimlist, int print_lvl = 0)
	{
		if (!empty && (OccPermGroup.rank+VirtPermGroup.rank+ppg.rank > 0))
		{
			dset = DsymToSets(dsym);
			ppg.clear(); 
			ppg = ProductPermutationGroup(rank, dset);
		}
		// VERY confusing but to decompose it... This function returns: 
		// a vector of (Cmap,(Operm,Vperm)) 
		std::vector<std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > > to_re;
		std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > to_push;		
		std::vector<int> blank; 
		for (int dim = 0; dim < rank/2 ; dim++)
			blank.push_back(dim);
		blank.push_back(1);
		if (empty || rank == 2 || (OccPermGroup.rank+VirtPermGroup.rank+ppg.rank == 0))
		{
			to_push = std::make_pair(dimlist,std::make_pair(blank,blank));
			to_re.push_back(to_push);
			return to_re; 
		}		
		//Construct 0-rank range which will be reused throughout. 
		std::vector<int> tmp, rang, tmp2;
		for (int dim = 0 ; dim < rank ; dim++)
			rang.push_back(dim);
		std::set< std::vector<int> > TmpSet;  // To reduce stupidity ;P 
		if (rank == 2 || OneGroup) //Dealt with in a specal way only one PG.
		{
			ppg.Begin();
			while(true)
			{
				tmp = rang; 
				MutVectorShiftInPlace(tmp , ppg.Current() ,0);	
				tmp2 = SuckDims(dimlist,tmp);
				if (TmpSet.count(tmp2) == 0)
				{
					to_push = std::make_pair(tmp2, std::make_pair(ppg.Current(),ppg.Current()));
					to_re.push_back(to_push);
					TmpSet.insert(tmp2); 
				}				
				if (!ppg.Next())
					break; 
			}			
			return to_re; 
		}
		if(OccPermGroup.size() == 0 || VirtPermGroup.size() == 0)
		{
			cout << "Can't Explode Cmap without OccPermGroup,VirtPermGroup " << OccPermGroup.size() << " :: " << VirtPermGroup.size() << endl;  
			system("sleep 10000");
		}
		to_re.reserve(OccPermGroup.size()*VirtPermGroup.size()); 
		//Explode Cdims... 
		OccPermGroup.Begin();
		while(true)
		{
			VirtPermGroup.Begin();
			while(true)
			{
				tmp = rang; 				
				MutVectorShiftInPlace(tmp,OccPermGroup.Current(),0);	
				MutVectorShiftInPlace(tmp,VirtPermGroup.Current(),rank/2);
				tmp2 = SuckDims(dimlist,tmp);
				//Figure out if this presents a new Contraction. 
				if (TmpSet.count(tmp2) == 0)
				{
					to_push = std::make_pair(tmp2, std::make_pair(OccPermGroup.Current(), VirtPermGroup.Current()));
					to_re.push_back(to_push);
					TmpSet.insert(tmp2); 
				}				
				if (!VirtPermGroup.Next())
					break; 
			}				
			if (!OccPermGroup.Next())
				break; 
		}		
		return to_re; 
	}

	// See routine above. 
	// Logic:  Loop pairs of dims in cdim1
	// are they in the same list of dset?  if so map to cdim2
	// are they in the same list of dset?  if so are these two in the same order. 
	// if not return false; 
	bool CdimsAreCompatible(const std::vector<int>& cd1, const std::vector<int>& cd2,
	const std::vector<std::vector<int> >& set1, const std::vector<std::vector<int> >& set2 )
	{
		if (!set1.size() || !set2.size())
			return true; 
		std::vector<std::vector<int> >::const_iterator iter; 
		std::vector<int>::const_iterator siter; 
		std::vector<int>::const_iterator vit1; 
		std::vector<int>::const_iterator vit2;
		bool found11; 
		int mapsto1 = 0 ; 
		int mapsto2 = 0 ;
		
		for (vit1 = cd1.begin(); vit1 + 1 != cd1.end(); ++vit1)
		{
			for (vit2 = vit1 + 1; vit2 != cd1.end(); ++vit2)
			{
				found11 = false; 
				for (iter = set1.begin(); iter != set1.end(); ++iter)
				{
					if (iter->size() == 0)
						continue; 
					siter = std::find(iter->begin(), iter->end(), *vit1); 
					if (siter != iter->end())
					{
						found11 = true; 
						break; 
					}
				}
				if (!found11)
					break; 
				if (std::find(iter->begin(), iter->end(), *vit2) == iter->end())
					break; 
				// Map onto cd2
				mapsto1 = cd2[(int)(vit1-cd1.begin())]; 
				mapsto2 = cd2[(int)(vit2-cd1.begin())]; 
				// Do these lie in the same set? 
				for (iter = set2.begin(); iter != set2.end(); ++iter)
				{
					if (std::find(iter->begin(), iter->end(), mapsto1) == iter->end())
						continue; 
					else 
					{
						if (std::find(iter->begin(), iter->end(), mapsto2) == iter->end())	
							break; 
						else
						{ 
							if ( ((*vit1 < *vit2) && (mapsto1 > mapsto2)) || ((*vit1 > *vit2) && (mapsto1 < mapsto2)) )
							{				 
								return false; 
							}
							else 
								break; 
						}
					}
				}
			}
		} 
		return true; 
	}

	// for the special case of Amplitude-Amplitude contractions. 
	bool CdimsAreCompatible(const std::vector<int>& cd1, const std::vector<int>& cd2, int rnk1, int rnk2, int print_lvl = 0)
	{
		std::vector<int> o1,o2;
		std::vector<int> v1,v2;
		int no1 = 0 ; 
		int no2 = 0 ;		
		for (std::vector<int>::const_iterator iter = cd1.begin(); iter != cd1.end(); ++iter)
		{
			if (*iter < rnk1/2)
				++no1; 
		}
		for (std::vector<int>::const_iterator iter = cd2.begin(); iter != cd2.end(); ++iter)
		{
			if (*iter < rnk2/2)
				++no2; 
		}
		o1.reserve(no1);
		v1.reserve(cd1.size()-no1);
		o2.reserve(no2);
		v2.reserve(cd2.size()-no2);
		for (std::vector<int>::const_iterator iter = cd1.begin(); iter != cd1.end(); ++iter)
		{
			if (*iter < rnk1/2)
				o1.push_back(*iter); 
			else 
				v1.push_back(*iter); 
		}
		for (std::vector<int>::const_iterator iter = cd2.begin(); iter != cd2.end(); ++iter)
		{
			if (*iter < rnk2/2)
				o2.push_back(*iter); 
			else 
				v2.push_back(*iter); 
		}
		
		if (print_lvl)
		{
			cout << "In CdimsAreCompatible, o1,o2,v1,v2 sizes: " << o1.size() << o2.size() << v1.size() << v2.size() << endl;  
		}
		
		if (o1.size())
		{
			int s1 = SignedSort(o1);
			int s2 = SignedSort(o2);
			if (s1 != s2)
				return false;
		}
		if (v1.size())
		{
			int s1 = SignedSort(v1);
			int s2 = SignedSort(v2);			
			if (s1 != s2)
				return false; 			
		}
		return true; 
	}
	
	// Given the contracted dimensions of 2,2 and the old eimap 
	// This constructs a new eimap with a new sign... 
	// If It's an Amplitude-Amplitude contraction this also
	// Takes care of the resulting factor. 
	inline void GetSignEimap(std::map<int,std::pair<int,int> >& eimap , std::pair<std::vector<int> , std::pair<std::vector<int>,std::vector<int> > >& cit1 
	, std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > >& cit2 , int& eff_sign , std::map<int,std::pair<int,int> >& eff_eimap, int print_lvl = 0)
	{
		//Useful Shorthand for rank of 1,2 
		int DimO1, DimO2;
		bool trivial1 = false; 
		bool trivial2 = false; 
		DimO1 = eff_eimap[0].first;
		DimO2 = eff_eimap[0].second;
		std::map<int,std::pair<int,int> >::iterator eiter; 
		
		trivial1 = (DimO1 == 2); //If the rank of either tensor involved is 2 then there is no permutational stuff to do. 
		trivial2 = (DimO2 == 2);
		
		bool onegrp1 = (DimO1 == ((cit1.second).first).size()-1 && DimO1 != 2);
		bool onegrp2 = (DimO2 == ((cit2.second).first).size()-1 && DimO2 != 2);
		
		if (print_lvl)
		{
			cout << "Entered GetSignEimap With Dims: " << DimO1 << " , " << DimO2 << " Input Map: " << endl; 
			cout << "Inducing Permutations of one,two: " << endl; 
			if (!trivial1)
			{
				printvec(((cit1.second).first));
				printvec(((cit1.second).second));
			}
			if (!trivial2)
			{		
				printvec(((cit2.second).first));
				printvec(((cit2.second).second));
			}
		}
		//Get ready.. 
		eff_eimap.clear(); 
		eff_sign = 1;
		if (!trivial1 && !onegrp1)
		{
			eff_sign *= ((cit1.second).first)[DimO1/2]; 
			eff_sign *= ((cit1.second).second)[DimO1/2];
		}
		else if (onegrp1 && !trivial1)
			eff_sign *= ((cit1.second).first)[((cit1.second).first).size()-1]; 
		
		if (!trivial2 && !onegrp2)
		{
			eff_sign *= ((cit2.second).first)[DimO2/2]; 
			eff_sign *= ((cit2.second).second)[DimO2/2];
		}
		else if (onegrp2 && !trivial2)
			eff_sign *= ((cit2.second).first)[((cit2.second).first).size()-1]; 
		
		int dimtemp; 
		//Permute each element of eimap according to this contraction's generating permutations. 
		for (eiter = eimap.begin(); eiter != eimap.end() ; ++eiter)
		{
			dimtemp = ((eiter->second).second);
			
			if ((eiter->second).first == 0)
			{
				if (!trivial1 && !onegrp1)
					MutSingleDim(dimtemp, cit1.second, DimO1);
				else if (onegrp1)
				{
					int dt = dimtemp; 
					dimtemp = (cit1.second).first[dt];
				}
			}
			if ((eiter->second).first == 1)
			{
				if (!trivial2 && !onegrp2)
					MutSingleDim(dimtemp, cit2.second, DimO2);
				else if (onegrp2)
				{
					int dt = dimtemp; 
					dimtemp = (cit2.second).first[dt];
				}
			}
			eff_eimap[eiter->first] = std::pair<int,int>((eiter->second).first, dimtemp);
		}	
		return; 
	}

	// Takes Eimap some external dimensions. Gives you (!!! complement of !!!) a list of 
	// dimensions on One which come from that vector
	void ExternalToArg(std::map<int,std::pair<int,int> >& ef_eimap , std::vector<int>& ExternalVector, std::vector<int>& tore, int OneOrTwo)
	{
		int DimO1 = tore[0];  
		std::vector<int>::const_iterator it; 
		std::list<int> tmp; 
		
		for (it = ExternalVector.begin(); it != ExternalVector.end(); ++it)
		{
			if (OneOrTwo == 1)
			{
				if (ef_eimap[*it].first == 0)
				{
					tmp.push_back(ef_eimap[*it].second); 
				}			
			}
			else if (OneOrTwo == 2)
			{
				if (ef_eimap[*it].first == 1)
				{
					tmp.push_back(ef_eimap[*it].second); 
				}			
			}
		}
		tore.clear(); 
		tore.reserve(DimO1-tmp.size()); 
		for (int di = 0 ; di < DimO1 ; ++di)
		{
			if (std::find(tmp.begin(),tmp.end(),di) == tmp.end())
				tore.push_back(di); 
		}
		return ; 
	}
	
	//Dim is the source dim which will be mapped onto it's permuted counterpart
	//Size is rank of tensor on which these permutations lie
	//Perms are permutations of O,V dimensions. 
	inline void MutSingleDim(int& dim, std::pair<std::vector<int>,std::vector<int> >& perms,int& size)
	{
		if (dim > size)
			cout << "dim > size in MutSingleDim (dim , size)" << dim << " , "<< size << endl; 
		if (size%2)
			cout << "Size%2 in MutSingleDim : " << size << endl; 
		
		if (dim < size/2)
		{
			dim = (perms.first)[dim];
		}
		else if (dim >= size/2)
		{
			int tempdim = dim - size/2; 
			dim = (perms.second)[tempdim] + size/2;		
		}
	}
	
	// Permute many dimensions of an vector (with shift to do this to V)
	// Does it by reference. 
	inline void MutVectorShiftInPlace(std::vector<int>& _arg,const std::vector<int>& mutation,int d_shift)
	{
		std::vector<int> to_re(_arg); 
		for(int dim=0 ; dim < d_shift ; dim++)
		{
			to_re[dim] = _arg[dim];
		}
		for(int dim=0; dim < (mutation.size()-1) ;dim++) //Last element of Mutation is the sign. 
		{
			to_re[dim+d_shift] = _arg[(mutation[dim])+d_shift];
		}
		_arg.swap(to_re); 
		return; 
	}
	
	inline long int fac(const int& number) const 
	{
		if (number > 12 || number < 0)
		{
			cout << "Please Add this case to TensorSymmetry.C" << endl;  
			return 0;  
		}
		switch (number)
		{
			case 0:
				return 1;
			case 1:
				return 1;
			case 2:
				return 2;
			case 3:
				return 6;
			case 4:
				return 24;
			case 5:
				return 120;
			case 6:
				return 720;
			case 7:
				return 5040;
			case 8:
				return 40320;
			case 9:
				return 362880;
			case 10:
				return 3628800;
			case 11:
				return 39916800;
			case 12:
				return 479001600;
		}
	}
	
	// Sort the Dimensions of a std::vector & return signature. 
	inline int SignedSort(std::vector<int>& array)  const 
	{
		if (array.size() == 2)
		{
			if (array[0] > array[1])
			{
				std::swap(array[0],array[1]); 
				return -1; 
			}
			else if (array[0]==array[1])
				return 0;
			else 
				return 1; 
		}
		int sign = 1; 
		std::vector<int>::iterator j;
		int flag = 1;    // set flag to 1 to begin initial pass 
		while (flag)
		{
			flag = 0;
			for(j = array.begin(); j+1 != array.end() ; ++j)
			{
				if (*(j+1) < *j)      
				{ 
					sign *= -1; 
					std::iter_swap(j,j+1); 
					flag = 1;              
				}
				else if (*(j+1) == *j)
					return 0; 
			}
		}
		return sign; 
	}
	
	// Sort the Dimensions of this index and return the sign 
	// Uses Bubblesort because the lists are short and the sign is easily obtained. 
	// Enforces a Pauli Antisymmetry. 
	// Fast logic for ranks 2,3. 
	inline int SortDex(std::vector<int>& array) const 
	{
		int Sz = array.size();
		int sign = 1; 
		if (Sz <= 2) 
			return 1; 
		else if (Sz == 4)
		{
			if (array[0] > array[1])
			{
				std::swap(array[0], array[1]);
				sign *= -1; 
			}
			else if (array[0] == array[1])
				return 0;
			if (array[2] > array[3])
			{
				std::swap(array[2], array[3]);
				sign *= -1; 			
			}
			else if (array[2] == array[3])
				return 0;
			return sign; 
		}
		int flag = 1;    // set flag to 1 to begin initial pass
		std::vector<int>::iterator i, j;
		i = array.begin() + Sz/2; 
		//Sort O 	
		while (flag)
		{
			flag = 0;
			for(j = array.begin(); j+1 != i ; ++j)
			{
				if (*(j+1) < *j)      
				{ 
					sign *= -1; 
					std::iter_swap(j,j+1); 
					flag = 1;              
				}
				else if (*(j+1) == *j)
					return 0; 
			}
		}
		//Sort V
		flag=1;
		while (flag)
		{
			flag = 0;
			for(j = i ; j+1 != array.end() ; ++j)
			{
				if (*(j+1) < *j)      
				{ 
					sign *= -1; 
					std::iter_swap(j,j+1); 
					flag = 1;              
				}
				else if (*(j+1) == *j)
					return 0; 	
			}
		}
		return sign; 
	}
	
	// Sort the Dimensions of this index 
	// But only do so according to permutations that are present in the index
	// As Perscribed by dset of this tensor. 
	inline int SortDex(std::vector<int>& array , std::vector<std::vector<int> >& syms) 
	{	
		if(syms.size() == 0 || rank == 2 || empty)
			return 1; 
		int sign = 1; 
		std::vector<std::vector<int> >::iterator Si;
		
		for (Si = syms.begin(); Si != syms.end(); ++Si)
		{
			if (Si->size() == 0)
				continue; 
			else
			{
				if (Si->size() == 2) // Optimal logic for this simple case of two Dimensions. 
				{
					if(array[(*Si)[0]] > array[(*Si)[1]])
					{
						std::swap(array[(*Si)[0]], array[(*Si)[1]]);    
						sign *= -1; 
					}
					 else if (array[(*Si)[0]] == array[(*Si)[1]])
						 return 0; 
				}
				else
				{
					std::vector<int> tvec(SuckDims(*Si,array));
					sign *= SignedSort(tvec);
					if (sign != 0)
						BlowDims(tvec,*Si,array); 
					else 
						return 0; 
				}
			}
		}
		return sign; 
	}

	// Sort the Dimensions of this index 
	// But only do so according to permutations that are present in the index
	// As Perscribed by dset of this tensor. 
	inline int SortDex(std::vector<int>& array , std::vector<std::vector<int> >& syms, const std::vector<int>& UnFacExto, const std::vector<int>& UnFacExtv) 
	{	
		if((syms.size() == 0 && !(UnFacExto.size() + UnFacExtv.size())) || rank == 2 || empty)
			return 1; 
		int sign = 1; 
		std::vector<std::vector<int> >::iterator Si;
		for (Si = syms.begin(); Si != syms.end(); ++Si)
		{
			if (Si->size() < 2)
				continue; 
			else
			{
				if (Si->size() == 2) // Optimal logic for this simple case of two Dimensions. 
				{
					if(array[(*Si)[0]] > array[(*Si)[1]])
					{
						std::swap(array[(*Si)[0]], array[(*Si)[1]]);
						sign *= -1; 
					}
					else if (array[(*Si)[0]] == array[(*Si)[1]])
						return 0; 
				}
				else
				{
					std::vector<int> tvec(SuckDims(*Si,array));
					sign *= SignedSort(tvec);
					if (sign != 0)
						BlowDims(tvec,*Si,array); 
					else 
						return 0; 
				}
			}
		}
		if(UnFacExto.size())
		{	
			if (UnFacExto.size() == 2) // Optimal logic for this simple case of two Dimensions. 
			{
				if (array[UnFacExto[0]] == array[UnFacExto[1]])
					return 0; 
				else if(array[UnFacExto[0]] > array[UnFacExto[1]])
				{
					std::swap(array[UnFacExto[0]], array[UnFacExto[1]]); 
					sign *= -1; 
				}
			}
			else if (UnFacExto.size() > 2)
			{
				std::vector<int> tvec(SuckDims(UnFacExto,array));
				sign *= SignedSort(tvec);
				if (sign != 0)
					BlowDims(tvec, UnFacExto, array); 
				else 
					return 0; 
			}
		}
		if (UnFacExtv.size())
		{
			if (UnFacExtv.size() == 2) // Optimal logic for this simple case of two Dimensions. 
			{
				if (array[UnFacExtv[0]] == array[UnFacExtv[1]])
					return 0; 
				else if(array[UnFacExtv[0]] > array[UnFacExtv[1]])
				{
					std::swap(array[UnFacExtv[0]], array[UnFacExtv[1]]); 
					sign *= -1; 
				}
			}
			else if (UnFacExtv.size() > 2)
			{
				std::vector<int> tvec(SuckDims(UnFacExtv,array));
				sign *= SignedSort(tvec);
				if (sign != 0)
					BlowDims(tvec, UnFacExtv, array); 
				else 
					return 0; 
			}
		}
		return sign; 
	}	
	
	// Ignores Negative one. 
	bool ObeysMyOrdering(const std::vector<int>& arg) const 
	{
		std::vector<int>::const_iterator viter,viter2; 
		std::vector<std::vector<int> >::const_iterator Si;		
		for (Si = dset.begin(); Si != dset.end(); ++Si)
		{
			if (Si->size() == 0)
				continue; 
			if (Si->size() == 2) // Optimal logic for this simple case of two Dimensions. 
			{
				if(arg[(*Si)[0]] < 0 || arg[(*Si)[1]] < 0)
					continue; 
				else if(arg[(*Si)[0]] > arg[(*Si)[1]])
					return false; 
			}
			else 
			{
				for (viter = Si->begin(); viter != Si->end(); ++viter)
				{
					for (viter2 = viter ; viter2 != Si->end(); ++viter2)
					{
						if (*viter == *viter2)
							continue;
						if (arg[*viter] < 0 || arg[*viter2] < 0)
							continue;
						else if (*viter < *viter2)
						{
							if (arg[*viter] > arg[*viter2])
								return false; 
						}
						else 
						{
							if (arg[*viter] < arg[*viter2])
								return false; 							
						}	
					}
				}
			}
		}
		return true; 
	}
	
	// Both Tensors's contracted dimensions are kept sorted in ExpandCpair(). 
	// This means the redundant pairs which only permute the order of 
	// contracted dimensions are omitted. 
	long int ObtainCdimFactor(std::vector<int>& ToSort, int print_lvl = 0 )
	{
		long int tore = 1;
		if (print_lvl)
		{
			cout << "Obtaining Cdim for: ";
			printvec(ToSort); 
		}
		std::vector<std::vector<int> >::const_iterator dsiter; 
		std::vector<int>::const_iterator viter; 
		for (dsiter = dset.begin(); dsiter != dset.end(); ++dsiter)
		{
			int NumberInThisSet = 0; 
			for (viter = dsiter->begin(); viter != dsiter->end(); ++viter)
			{
				if (ToSort[*viter] >= 0)
					++NumberInThisSet; 
			}
			if (print_lvl)
			{
				cout << "Num in this set: " << NumberInThisSet << " : ";
				printvec(*dsiter); 
			}
			long int Catch1 = tore; 
			long int Catch2 = fac(NumberInThisSet); 
			long int Catch3 = INT_MAX; 
			if (Catch1*Catch2 > Catch3)
			{
				cout << " !*!*!*!*!*!*!*!*!*!*!*!  Danger of integer rollover. !*!*!*!*!*!*!*!*!*!*!*!*!" << endl;  
				cout << Catch1*Catch2 << " " << Catch3 << endl; 
			}
			tore *= fac(NumberInThisSet); 
		}
		return tore; 
	}

	// I apologize in advance for the existence of this routine. 
	// Read aloud slowly. And Repeat. You'll figure it out. 
	// -1 -1 1 0 => 3 2
	std::vector<int> OrderOfDimensionsToDimensionsInOrder(const std::vector<int>& arg) 
	{	
		std::list<int> tmp; 
		std::vector<int>::const_iterator pI2; 
		for (int p = 0 ; p < arg.size() ; ++p)
		{
			pI2 = std::find(arg.begin(), arg.end(), p); 
			if (pI2 == arg.end())
				break; 
			else 
				tmp.push_back(pI2 - arg.begin());
		}
		std::vector<int> tore(tmp.size());
		tore.assign(tmp.begin(), tmp.end());
		return tore; 
	}
	
	// Are the dimensions of this vector sorted (in ascending order) 
	bool IsPlainlySorted(const std::vector<int>& arg) const 
	{
		std::vector<int>::const_iterator i, j ;
		for (i = arg.begin(); i != arg.end(); ++i)
		{
			for (j = i+1; j != arg.end(); ++j)
			{
				if (*i < 0 || *j < 0)
					continue; 
				if (*i > *j)
					return false; 
			}		
		}
		return true; 
	}
	
	// Are the dimensions of this index sorted (in ascending order) 
	// Is this Called at all? 
	bool IsSorted(const std::vector<int>& arg) const 
	{
		if(rank == 2)
			return true; 
		std::vector<int>::const_iterator i, j ;
		i = arg.begin() + arg.size()/2; 
		for(j = arg.begin(); j+1 != i ; ++j)
		{
			if (*(j+1) <= *j)      
				return 0; 
		}
		//Sort V
		for(j = i ; j+1 != arg.end() ; ++j)
		{
			if (*(j+1) <= *j)      
				return 0; 	
		}
		return true; 
	}
	
	// Are the dimensions of this index sorted (in ascending order) ? 
	// According to the reduced symmetry as specified by syms? 
	bool IsSorted(const std::vector<int>& arg ,const std::vector<std::vector<int> >& syms) const 
	{
		if(rank == 2 || empty)
			return true; 
		std::vector<std::vector<int> >::const_iterator it1;
		std::vector<int>::const_iterator it2,it3;
		for (it1 = syms.begin() ; it1 != syms.end() ; ++it1)
		{
			for(it2 = it1->begin() ; it2 != it1->end() - 1 ; ++it2)
			{
				it3 = it2; 
				it3++; 
				if (it3 == it1->end())
					break; 
				if (arg[*it2] >= arg[*it3])
					return false; 
			}
		}
		return true; 
	}
	
	inline void ShiftVec(std::vector<int>& arg,int shift)
	{
		std::vector<int>::iterator ait; 
		for (ait = arg.begin() ; ait != arg.end() ; ait++)
			(*ait) += shift; 
		return; 
	}

	// Order of Product Permutation Group Specified by Dsym. 
	int GroupSizeForDsym(std::map<int,std::vector<int> >& Dsym) const
	{
		int to_re = 1; 
		std::vector<std::vector<int> > Tmp = DsymToSets(Dsym);
		for (std::vector<std::vector<int> >::const_iterator t = Tmp.begin() ; t != Tmp.end(); ++t) 
			to_re *= fac(t->size()); 
		return to_re;  
	}
	
	// Collect the unique subsets of symmetrical dimensions in dsym 
	std::vector<std::vector<int> > DsymToSets(std::map<int,std::vector<int> >& Dsym) const
	{
		std::vector<int> tmp; 
		std::vector<std::vector<int> > to_re; 
		to_re.reserve(rank);
		tmp.reserve(rank);
		std::map<int,std::vector<int> >::iterator si;
		for (si = Dsym.begin(); si != Dsym.end() ; ++si)
		{
			tmp = si->second; 
			if (tmp.size() == 0)
				continue; 
			tmp.push_back(si->first); 
			std::sort(tmp.begin(), tmp.end());
			Unique(tmp);
			to_re.push_back(tmp); 
		}
		Unique<std::vector<std::vector<int> > >(to_re);
		// Remove any empty sets; 
		for (std::vector<std::vector<int> >::iterator ti = to_re.begin(); ti != to_re.end(); ++ti)
		{	
			if (ti->size() == 0)
				to_re.erase(ti);
		}
		return to_re;  
	}
	
	// Specific for first-half last-half type symmetry. 
	inline void DsymToDimDim(std::map<int,std::vector<int> >& Dsym, std::vector<std::pair<int,int> >& cosym, std::vector<std::pair<int,int> >& cvsym)
	{
		cosym.clear();
		cvsym.clear(); 
		cosym.reserve(rank); 
		cvsym.reserve(rank);
		std::vector<std::vector<int> > Sets = DsymToSets(Dsym);
		std::vector<std::vector<int> >::iterator s0; 
		std::vector<int>::iterator s1,s2; 
		for (s0 = Sets.begin() ; s0 != Sets.end() ; ++s0)
		{
			for ( s1 = s0->begin() ; s1 != s0->end() ; ++s1)
			{
				for (s2 = s1 ; s2 != s0->end() ; ++s2)
				{
					if (*s1 < *s2)
					{
						if(*s1 < rank/2)
						{
							cosym.push_back(std::make_pair(*s1,*s2));
						}
						else if (*s1 >= rank/2)
						{
							cvsym.push_back(std::make_pair(*s1-rank/2,*s2-rank/2));
						}
					}
				}
			}	
		}
	}

	// Are the two indices in order in the arguement permutation. 
	inline bool Mut(std::vector<int>& arg,int first, int second)
	{
		int fplace=-1;
		int splace=-1;
		int tmp; 
		for(int dim = 0 ; dim < (arg.size()-1) ; dim++)
		{
			if (arg[dim] == first)
				fplace=dim;
			else if (arg[dim] == second)
				splace=dim;
		}
		return (fplace < splace); 
	}
	
	void printvec(const std::vector<int>& arg, bool NoLineBreak = false) const
	{	
		std::vector<int>::const_iterator pI2; 
		for (pI2 = arg.begin() ; pI2 != arg.end() ; ++pI2)
		{
			cout << " " << *pI2; 
		}
		if (!NoLineBreak)
			cout << endl; 
	}
	
	void printvov(const std::vector<std::vector<int> > arg) 
	{	
		std::vector<std::vector<int> >::const_iterator pI1; 
		std::vector<int>::const_iterator pI2; 
		std::vector<int> temp; 
		
		for (pI1 = arg.begin() ; pI1 != arg.end() ; ++pI1)
		{
			temp = *pI1; 
			for (pI2 = temp.begin() ; pI2 != temp.end() ; ++pI2)
			{
				cout << " " << *pI2;
			}
			cout << endl; 
		}
		cout << endl; 
	}
	
	// Debug routine. Display Cdims and Generating Permutations. 
	void PrintCdims(std::vector<std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > >& arg) const
	{
		cout << "Cdim Sz: " << arg.size() << endl; 
		std::vector<std::pair<std::vector<int>,std::pair<std::vector<int>,std::vector<int> > > >::const_iterator iter; 
		for (iter = arg.begin(); iter != arg.end(); ++iter)
		{
			printvec(iter->first, true);
			cout << " || " ; 
			printvec(iter->second.first, true); 
			cout << " : " ; 
			printvec(iter->second.second, true);
			cout << endl; 
		}
		return;
	}
	
	void Print()
	{
		cout << " -- Symmetry Information -- " << endl;
		cout << "Rank: " << rank;
		cout << " Dsym: " << dsym.size() << endl;
		std::map<int, std::vector<int> >::iterator iter; 
		for (iter = dsym.begin() ; iter != dsym.end(); ++iter)
		{
			cout << iter->first << " :: " ;	
			printvec(iter->second); 
		}
		cout << "O group : " <<  endl; 
		OccPermGroup.Print();
		cout << "V group : " <<  endl; 
		VirtPermGroup.Print(); 
		cout << "PPG : " <<  endl; 
		ppg.Print(); 
	}
	
};

#endif
#endif 
