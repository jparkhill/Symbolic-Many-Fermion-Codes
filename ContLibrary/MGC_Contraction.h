#ifndef NOMGC

#include <string>
#include <vector>
#include <map>
#include "blas_include.h"

#include "MGC_Dense.h"
#include "MGC_Timer.h"
#include "MGC_Sets.h"
#include "MGC_MBTensor.h"
#include "MGC_filtering.h"
#include "MGC_StoreRecall.h"

#ifndef _contraction_h
#define _contraction_h

//extern TimeProfile Profiler;
//extern StatCollector Stats;
using namespace std; 

////////////////////////////////////////////////////////////////////////////////
// contraction, tensor products etc. 
// They suppose the functions on a tensor from MGC_GenTensor.C
// and a TensorSymmetry object.  
// 
// To use Contract(); just populate it, then execute it. cf. MGC_ampeqns.C
//
// John Parkhill 2008;
////////////////////////////////////////////////////////////////////////////////

// To avoid recalculation of ExpandedContractionMap 
// This ContractionMapSource provides the recipe for each result. 
// so the the ExpandedContractionMap may be indexed. 
class ContractionMapSource
{
public:
	TensorSymmetry S1,S2; 
	std::vector<int> Cd1,Cd2; 
	ContractionMapSource(const TensorSymmetry A,const TensorSymmetry B,const std::vector<int>& C,const std::vector<int>& D)
	{
		S1.CopyFrom(A);	S2.CopyFrom(B); Cd1 = C; Cd2 = D; 
	}
	bool operator==(const ContractionMapSource& other)
	{
		return ((S1==other.S1)&&(S2==other.S2)&&(Cd1==other.Cd1)&&(Cd2==other.Cd2));
	}
};

// Simple container for the symmetrically unique version of a contraction
class ExpandedContractionMap : public std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >  
{
public:
	ExpandedContractionMap()
	{}
	void CopyFrom(const ExpandedContractionMap& other)
	{
		*this = other; 
	}
	ExpandedContractionMap& operator=(const std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >& other)
	{
		clear();
		resize(other.size());
		std::copy(other.begin(),other.end(),begin());
		return *this; 
	}
	void Write(Streamer& AStream)
	{
		if ( !(AStream.os->good()) )
		{
			cout << "Can't read from bad stream asfadgwwe" << endl;  
			system("sleep 10");
		}
		*(AStream.os) << " ExContMap { ";
		for (ExpandedContractionMap::const_iterator iter = begin(); iter != end(); ++iter)
		{
			AStream.WriteToStream(iter->first); 
			AStream.WriteToStream(iter->second.first); 
			AStream.WriteToStream(iter->second.second);
			AStream.Flush(); 
		}
		*(AStream.os) << " } ";
	}
	void Read(Streamer& AStream)
	{
		clear(); 
		std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > to_push; 
		if ( !(AStream.is)->good() )
			cout << "Can't read from bad stream asfsfee" << endl;  
	// read this block. 
	//	int strt = AStream.is->tellg();
	//	AStream.is->seekg (0, ios::end);
	//	int len = AStream.is->tellg();
	//	AStream.is->seekg (strt, ios::beg);
	//	char Buf[len];
	//	AStream.is->read(Buf,len);
	//		cout << "ReadBufLength:" << len << endl;  
		string sstart("ExContMap");
		string Tmp = AStream.ReadNextString(); 
		if (sstart != Tmp)
		{
			cout << "Does not begin Excontmap" << endl;
			cout << Tmp; 
		}
		// Read the vectors. 
		string bufs = AStream.ReadNextBraKetPair(); 
		int len = bufs.size(); 
		char Buf[len];
		std::copy(bufs.begin(), bufs.end() ,Buf);		
		for (int k = 0; k < len-5 ;)
		{
			to_push.first=ReadVector(Buf,k,len);
			to_push.second.first=ReadVector(Buf,k,len);
			to_push.second.second=ReadVector(Buf,k,len);
			push_back(to_push);
		}
	}
	char NextNonWhitespace(char* Buf, int i, int len)
	{
		for(int k = i ; k < len ; ++k)
		{
			if (Buf[k] > 32)
				return Buf[k];
		}
		return -1;
	}
	std::vector<int> ReadVector(char* Buf, int& i, int len)
	{
		std::vector<int> tore; 
		//		cout << "Reading vec from: " << len <<  " ";  
		//		cout.write(Buf+i,len);
		//		cout << endl << "------------------" << endl;  
		Buf += i ; 
		std::pair<int,int> bnds = NextBraketPair(Buf,i,len); 
		//	cout << endl << "bds: " << bnds.first << " " << bnds.second << endl;  
		//	cout.write(Buf+bnds.first,bnds.second-bnds.first);
		//	cout << endl; 
		std::list<std::string> tocast; 
		std::string tmp; int tmpi; 
		int k=bnds.first; 
		for (; k < bnds.second ; ++k) 
		{
			tmpi = int(Buf[k]); 
			//	cout << tmpi << ":" ;
			if (tmpi == 125)
			{
				if (ContainsASCIIint(tmp))
				{
					tocast.push_back(tmp);
					tmp.clear(); 
				}
				break;
			}
			if ((48 <= tmpi && tmpi <= 57) || tmpi == 45)
			{
				tmp.push_back(Buf[k]); 
			}
			else 
			{
				if (ContainsASCIIint(tmp))
				{
					tocast.push_back(tmp);
					tmp.clear(); 
				}
				else 
				{
					tmp.clear(); 
				}
			}
		}
		if (ContainsASCIIint(tmp))
			tocast.push_back(tmp);
		//	cout << "List of "  << tocast.size() << endl; 
		tore.reserve(tocast.size());
		for (std::list<std::string>::iterator iter = tocast.begin(); iter != tocast.end(); ++iter)
			tore.push_back(castint(*iter));
		i+=k+1; 
		//		cout << "ReadVector: " ;
		//		printvec(tore); 
		return tore; 
	}
	bool ContainsASCIIint(std::string& arg)
	{
		if (!arg.length())
			return false; 
		//	cout << "Contains: " << arg << endl;  
		int c=0; 
		for (std::string::iterator iter = arg.begin(); iter != arg.end(); ++iter)
		{
			c = int(*iter);
			if (48 <= c && c <= 57)
				return true; 
		}
		return false; 
	}
	int castint(string& arg)
	{
		//	cout << "Casting String: " << arg << endl;  
		int tore=0; 
		int c=0; 
		for (int dd=0 ; dd<arg.length() ; ++dd)
		{
			c=arg[arg.length()-(dd+1)];
			if (48 < c && c < 58)
				tore += (c-48) * pow(10.0,dd);
			else if (c == 45)
				tore *= -1; 
		}
		return tore; 
	}
	// exclusive. 
	std::pair<int,int> NextBraketPair(char* Buf, int i, int len)
	{
		//	cout << "Bra Ket strt" << i << " " << len << endl;  
		std::pair<int,int> tore;
		tore.second = len;
		for (int k = 0; k < len; ++k)
		{
			if (int(Buf[k] == 123))
				tore.first = k; 
			else if (int(Buf[k] == 125))
			{
				tore.second = k; 
				break; 
			}
		}
		return tore; 
	}
	void print()
	{
		std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >::const_iterator MWrite; 
		cout << "Cmap: " << endl;  
		for (MWrite = begin(); MWrite != end(); ++MWrite)
		{
			printvec(MWrite->first,true); 
			cout << " : ";
			printvec(MWrite->second.first,true); 
			cout << " , " ; 
			printvec(MWrite->second.second,true);
			cout << endl;   
		}
	}	
	inline void printvec(const std::vector<int> arg, bool IfNoLine = false)
	{	
		std::vector<int>::const_iterator pI2; 
		for (pI2 = arg.begin() ; pI2 != arg.end() ; ++pI2)
		{
			cout << " " << *pI2; 
		}
		if (!IfNoLine)
			cout << endl; 
	}
};


// Contains parameters and things directly derived from parameters of 
// a contraction by things like symmetry information. 
// These are packaged and kept around to avoid the overhead of recomputation 
// per iteration for > rank 4 models. 
class Contraction
	{
	public: 
		bool populated; // Basic Parameters inserted		
		bool prepared;  // Derived parameters initalized 
		bool executed;  // contracted. 
		//Contraction Paramters
		std::map<int,pair<int,int> > eimap;
		std::map<int,pair<int,int> > cmap;
		double factor;
		int algorithm;  // Multipliers require alg = 0. Alters if intermediates are given complete bra-ket sym. 
						// Algorithm = 3 alters the way blocking is done. 
		std::vector<int> external_eindices; 
		//Derived Quantities
		std::map<int, std::vector<int> > osym; //Symmetries of output.. 
		std::vector<std::vector<int> > oset; 
		std::vector<int> Cd1,Cd2; 
		ExpandedContractionMap cmap1,cmap2;   
		std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > Cpair1, Cpair2; 
		std::vector<int> o_type;// Orbital type of the result.
		std::vector<int> BraketType; // 1 Up // -1 Down 
		TensorSymmetry SymmetryOfOutput; 
		Filters FilterOfOutput;
		
		Contraction() : factor(1.0), algorithm(0), populated(false), prepared(false), executed(false) 
		{ return; }
		Contraction(const Contraction& other) : eimap(other.eimap) , cmap(other.cmap), BraketType(other.BraketType),
			factor(other.factor), algorithm(other.algorithm), external_eindices(other.external_eindices),
			Cd1(other.Cd1), Cd2(other.Cd2), cmap1(other.cmap1), cmap2(other.cmap2), o_type(other.o_type), SymmetryOfOutput(other.SymmetryOfOutput), 
			FilterOfOutput(other.FilterOfOutput), populated(other.populated), prepared(other.prepared), executed(other.executed) 
		{
			return; 
		}
		void Write(Streamer& Dest)
		{
			cmap1.Write(Dest); cmap2.Write(Dest); Dest.Flush();  
		}
		void Read(Streamer& Dest)
		{
			cmap1.Read(Dest); cmap2.Read(Dest);  
		}
		void Populate(std::map<int,pair<int,int> >& a_eimap , std::map<int,pair<int,int> >& a_cmap , double a_factor , 
					  int a_algorithm , std::vector<int> a_external_eindices)
		{
			if (populated || prepared || executed)
				return;
			eimap = a_eimap;
			if (a_cmap.size())
				cmap = a_cmap;
			o_type = std::vector<int>(a_eimap.size(),0);
			factor = a_factor;
			algorithm = a_algorithm; 
			external_eindices = a_external_eindices; 
			populated = true; 
		} 
		// Unused? JAP 12/09
		void Populate(std::map<int,pair<int,int> >& a_eimap , std::map<int,pair<int,int> >& a_cmap, bool a_if_anti , double a_factor , 
					  int a_algorithm , std::vector<int> a_external_eindices)
		{
			if (populated || prepared || executed)
				return;
			eimap = a_eimap;
			if (a_cmap.size())
				cmap = a_cmap;
			o_type = std::vector<int>(a_eimap.size(),0);
			factor = a_factor;
			algorithm = a_algorithm; 
			external_eindices = a_external_eindices; 
			populated = true; 
		}

		// Translates arguement tensors into the output. 
		// Most costly step is: ExpandCpair/ExplodeCPair
		// which can be avoided by Read(); 
		//
		// If the intermediate shares > 2 lines with another amplitude
		// Then enforcing Bra-Ket symmetry on the intermediate (As Occurs if SortInter = true)
		// will result in an adjustment for 
		// The factor of the output which can not yet be determined in this code. 
		// This is controlled by IfSortInter.  (JAP 12/2009)
		void PrepareContraction(GenTensor& one, GenTensor& two, GenTensor& three, bool SortInter = false, int print_lvl = 0)
		{
			SortInter=false; 
			// although cmap can be re-used these must be transmitted to the
			// intermediate. 
			if (prepared && three.IsInter)
			{	
				three.MySym.CopyFrom(SymmetryOfOutput);
				three.BraketType = BraketType;
				three.o_type = o_type; 
				three.num_pairs = one.num_pairs; 
				// Transmit for blocking. 
				three.MyTypes = one.MyTypes; 
				if (oset.size())
				{
					three.MySym.rank = three.Rank(); 
					three.MySym.AssignPG(osym);
				}
				else 
				{				
					three.MySym.empty = true; 
					three.MySym.rank = three.Rank(); 
				}	
				return; 
			}
			else if (prepared)
				return; 
			// Otherwise begin anew. 
			osym.clear(); oset.clear();Cd1.clear();Cd2.clear(); 
			// cmap1.clear();cmap2.clear();
			o_type = std::vector<int>(eimap.size(),0);
			SymmetryOfOutput.clear(); 
			BraketType = std::vector<int>(eimap.size(),0);
			// --------------- Begin Symmetry mapping section -----------------
			map<int,pair<int,int> >::iterator ei_iter; 	
			std::vector<int> tempv; 
			std::vector<int>::iterator tempvi;
			std::map<int,int> o_dims_in_o,t_dims_in_o; // Maps (dim of one) -> dim of output
			std::map<int,int>::iterator odimi; 
			if (three.Rank() != 2)
			{
				//Map these onto the symmetries of output. 
				for (ei_iter = eimap.begin() ; ei_iter != eimap.end() ; ei_iter++)
				{
					if((ei_iter->second).first == 0) 
					{
						o_dims_in_o[((ei_iter->second).second)] = ei_iter->first; 
					}
					else 
						t_dims_in_o[((ei_iter->second).second)] = ei_iter->first; 
				}
				for(int rdim = 0 ; rdim < three.Rank() ; rdim++)
				{
					if(eimap[rdim].first == 0) //Symmetries come from one. 
					{
						if (three.IsInter) //Write output orbital type. 
						{
							three.o_type[rdim] = one.o_type[(eimap[rdim].second)];
							o_type[rdim] = one.o_type[(eimap[rdim].second)];
						}
						if (one.Rank() == 2)
						{
							osym[rdim];
							continue; 
						}
						tempv = one.MySym.dsym[(eimap[rdim].second)];
						for (tempvi = tempv.begin() ; tempvi != tempv.end() ; tempvi++) //look to see if each index this is sym to is in output
						{
							odimi = o_dims_in_o.find(*tempvi);
							if (odimi != o_dims_in_o.end())
								osym[rdim].push_back(o_dims_in_o[*tempvi]);
						} 
					}
					else //Symmetries come from two.
					{
						if (three.IsInter) //Write output orbital type. 
						{
							three.o_type[rdim] = two.o_type[(eimap[rdim].second)];
							o_type[rdim] = two.o_type[(eimap[rdim].second)];
						}
						if (two.Rank() == 2)
						{
							osym[rdim];
							continue; 
						}
						tempv = two.MySym.dsym[(eimap[rdim].second)];
						for (tempvi = tempv.begin() ; tempvi != tempv.end() ; tempvi++) //look to see if each index this is sym to is in output
						{
							odimi = t_dims_in_o.find(*tempvi);
							if (odimi != t_dims_in_o.end())
								osym[rdim].push_back(t_dims_in_o[*tempvi]);
						} 		
					}
				}
			}
			if (three.IsInter) //Compute output orbital type(s) & Symmetry & Filter for an intermediate
			{
				// Differentiate Bra-Ket indices for the purposes of symmetry
				for(int rdim = 0 ; rdim < three.Rank() ; rdim++)
				{
					if(eimap[rdim].first == 0) //Symmetries come from one. 
					{
						three.o_type[rdim] = one.o_type[(eimap[rdim].second)];
						o_type[rdim] = one.o_type[(eimap[rdim].second)];
						three.BraketType[rdim] = one.BraketType[(eimap[rdim].second)];
						BraketType[rdim] = one.BraketType[(eimap[rdim].second)];						
					}
					else
					{
						three.o_type[rdim] = two.o_type[(eimap[rdim].second)];
						o_type[rdim] = two.o_type[(eimap[rdim].second)];
						three.BraketType[rdim] = two.BraketType[(eimap[rdim].second)];
						BraketType[rdim] = two.BraketType[(eimap[rdim].second)];						
					}
				}			
				// Transmit for blocking.
				three.num_pairs = one.num_pairs; 
				three.MyTypes = one.MyTypes; 
				// Assign complete intermediate symmetry? 
				if (SortInter)
				{
					std::vector<int> ODownList; 
					std::vector<int> VDownList; 
					std::vector<int> OUpList; 
					std::vector<int> VUpList; 
					for(int rdim = 0 ; rdim < three.Rank() ; rdim++)
					{
						if (o_type[rdim] == 0)
						{
							if (BraketType[rdim] == -1)
								ODownList.push_back(rdim); 
							else 
								OUpList.push_back(rdim);
						}
						else 
						{
							if (BraketType[rdim] == -1)
								VDownList.push_back(rdim); 
							else 
								VUpList.push_back(rdim);
						}
					}
					// Bra indices not included. 
					if (ODownList.size() + OUpList.size() + VDownList.size() + VUpList.size())
					{
						for(int rdim = 0 ; rdim < three.Rank() ; rdim++)
						{
							if (o_type[rdim] == 0)
							{	
								if (BraketType[rdim] == -1)
									osym[rdim] = ODownList;
								else
									osym[rdim] = OUpList;
							}
							else 
							{	
								if (BraketType[rdim] == 1)
									osym[rdim] = VUpList;
								else
									osym[rdim] = VDownList;
							}
						}
					}
					else 
						osym.clear(); 
				}
				oset = three.MySym.DsymToSets(osym);
				if (oset.size())
				{
					// TO DO : COLLAPSE WITH Three.MySym. 
					SymmetryOfOutput.rank = three.Rank(); 
					three.MySym.rank = three.Rank(); 
					SymmetryOfOutput.AssignPG(osym);
					three.MySym.AssignPG(osym);
				}
				else 
				{
					for (int dim = 0 ; dim < three.Rank() ; dim++)
						osym[dim]=std::vector<int>();			
					SymmetryOfOutput.empty = true; 					
					three.MySym.empty = true; 
				}				
				if (print_lvl)
				{
					cout << " -------------SymmetryOfOutput(PCDBG)------------ " << endl; 
					for (std::vector<std::vector<int> >::iterator tm = oset.begin(); tm !=oset.end() ;++tm)
					{
						for (std::vector<int>::iterator t = tm->begin(); t != tm->end() ; ++t)
						{
							cout << " " << *t; 
						}							
						cout << endl; 
					}	
					SymmetryOfOutput.Print(); 
				}
			}
			else if (three.IsAmp)
			{
				o_type = three.o_type; 
				FilterOfOutput.push_back(new SpinFilter(three.Rank() , three.o_type , three.MyTypes.ospin, three.MyTypes.vspin)); 
			}
			else if (three.IsIntegral)
				SymmetryOfOutput.CopyFrom(three.MySym); 
			// TODO : External Sort for density Matrices. 
			
			// --------------- End Symmetry mapping section -----------------
			
			Cd1.clear(); 
			Cd2.clear(); 
			std::map<int,pair<int,int> >::const_iterator cmiter;			
			//Convert Cmap -> Dims in 1, dims in 2. 
			for(cmiter = cmap.begin(); cmiter != cmap.end(); ++cmiter)
			{
				Cd1.push_back(cmiter->second.first); 
				Cd2.push_back(cmiter->second.second); 
			}			
			
			if (!cmap1.size() && !cmap2.size())
			{
				if (Cd1.size() && one.IsAmp && (one.Rank() > 2) )
					one.MySym.ExpandCpair(Cd1, Cd2, two.MySym, cmap1, cmap2, print_lvl);
				else if (Cd1.size() && two.IsAmp && (two.Rank() > 2) )
					two.MySym.ExpandCpair(Cd2, Cd1, one.MySym, cmap2, cmap1, print_lvl);
				else 
				{ 
					if (print_lvl)
						cout << " .... Exploding CDims .... " << endl; 
					cmap1 = one.MySym.ExplodeCdims(Cd1, print_lvl);  //Cmap --> vector of all Possible Contractions in One and Corresponding Permutations
					cmap2 = two.MySym.ExplodeCdims(Cd2, print_lvl);  
				}

			}
			
			if (print_lvl>1)
			{
				cout << "Package Populated... " << endl; 
				cout << "Osym size: " << osym.size() << endl; 
				cout << "Oset size: " << oset.size() << endl; 
				cout << "Cd1 size: " << Cd1.size() << endl; 		
				cout << "Cd2 size: " << Cd2.size() << endl; 		
				cout << "Cmap1 size: " << cmap1.size() << endl; 		
				cout << "Cmap2 size: " << cmap2.size() << endl; 			
				cout << "Output Symmetries: " << endl; 
				SymmetryOfOutput.Print(); 
				cout << "Output Filter Length: " << FilterOfOutput.size() << endl; 
			}
			prepared = true; 
			return; 
		}
		
		// Algorithm determines nature of symmetry on intermediates. (c.f. above)
		void Contract(GenTensor& one, GenTensor& two, GenTensor& three, int print_lvl = 0, std::vector<int> dbgdex = std::vector<int>())
		{
			if (print_lvl>1)
				printf("Pop = %d, Exe = %d, Alg = %d\n", populated , executed, algorithm);
			if (!populated)
			{
				cout << "Error!!! Contracting Unpopulated package" << endl; 	
				system("sleep 100"); 
			}

			// 3 is for NoPairBlocking
			if (algorithm == 0)
				PrepareContraction(one, two, three, false, print_lvl);
			else if (algorithm == 1 || algorithm == 3)
				PrepareContraction(one, two, three, true, print_lvl);
			else 
				cout << "Invalid Algorithm..  " << endl; 
			
			if (three.IsIntegral)
				ContractLayeredDM(one, two, three); 
			else if(three.IsAmp)
				ContractLayered(one, two, three, print_lvl); 				
			else if (three.IsInter) 
			{
				ContractLayeredInter(one, two, three, print_lvl);		
				three.ReStruc(pow(10.0,-15.0));
			}				 
			executed = true; 
		}
		
		void ContractLayered(GenTensor& one, GenTensor& two, GenTensor& three , int print_lvl = 0 , std::vector<int> DBGDEX = std::vector<int>())
		{
			if (!one.l2norm() || !two.l2norm())
				return;
			std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >::iterator cit1;
			std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >::iterator cit2;
			map<int,pair<int,int> > eff_eimap; //Effective Eimap. 
			int eff_sign=1;  //Effective sign
			
			DenseLayeredRep* dcr1;
			DenseLayeredRep* dcr2;
			
			if (algorithm != 3)
			{
				dcr1 = one.GetDLR();
				dcr2 = two.GetDLR(one.MyTypes);
			}
			else 
			{
				char* Tmp = "NoPairBlocking"; 
				std::vector<char*> Tmp2(1);
				Tmp2[0] = Tmp;
				dcr1 = one.GetDLR(Tmp2);
				dcr2 = two.GetDLR(Tmp2);
			}
			
			if(dcr1 == dcr2)
			{
				cout << "ERROR ONE==TWO in ContractLayered" << endl; 
				system("sleep 100");
			}
			
			double* dat1(NULL); 
			double* dat2(NULL); 
			int* Col1(NULL);
			int* Row1(NULL);
			int* Col2(NULL);
			int* Row2(NULL);
			dcr1->AllocForPartition((cmap1.begin())->first , Col1, Row1, dat1); 
			dcr2->AllocForPartition((cmap2.begin())->first , Col2, Row2, dat2); 
			
			//Tries to reduce function evalutation. 
			int Rank1 = one.Rank(); 
			int Rank2 = two.Rank(); 
			int Width1 = cmap1.begin()->first.size();
			int Width2 = cmap2.begin()->first.size(); 
			long int ColShift1 = 0; 
			long int ColShift2 = 0; 
			long int RowShift1 = 0; 
			long int RowShift2 = 0; 
			long int BlockBatchStart1 = 0; 
			long int BlockBatchStart2 = 0;
			long int BlockBatchSize1 = 0; 
			long int BlockBatchSize2 = 0; 

			if (print_lvl>2)
			{
				cout << "ContractLayered!" << endl; 
				one.MySym.PrintCdims(cmap1); 
				one.MySym.PrintCdims(cmap2); 
			}

			if (dcr1->clength*cmap2.size() > dcr2->clength*cmap1.size())
			{
				if (print_lvl)
					cout << "Len1 > Len2, picking first branch" << endl;   

				//Loop over possible contractions owing to permutational symmetry. 
				for (cit1 = cmap1.begin() ; cit1 != cmap1.end() ; ++cit1)
				{
					dcr1->ReBlock(cit1->first, true);
					if (print_lvl>1)
					{
						cout << "---------Blocking of 1-----------" << endl;  
						dcr1->MyBlocking.print(1);
					}
					for (cit2 = cmap2.begin() ; cit2 != cmap2.end() ; ++cit2)
					{
						if (!one.MySym.empty && !two.MySym.empty)
						{
							if (cit2 != cmap2.begin())
							{
								if (!one.MySym.OneGroup && !two.MySym.OneGroup)
								{
									if (!one.MySym.CdimsAreCompatible(cit1->first, cit2->first, one.Rank(), two.Rank()))
										continue; 
								}		
								else if (!one.MySym.CdimsAreCompatible(cit1->first, cit2->first, one.MySym.dset, two.MySym.dset))
									continue; 
							}
						}
						dcr2->ReBlock(cit2->first, false, print_lvl);
						if (print_lvl>1)
						{
							cout << "---------Blocking of 2-----------" << endl;  
							dcr2->MyBlocking.print(1);	
						}
						//Construct sign,eimap for this pair of contractions
						//Cheap way to pass Dim1 and Dim2
						eff_eimap.clear(); 
						eff_eimap[0] = std::pair<int,int>(one.Rank() , two.Rank() );
						three.MySym.GetSignEimap(eimap , *cit1 , *cit2 , eff_sign , eff_eimap , print_lvl);
						
						dcr1->MyBlocking.Begin();
						while (true)
						{
							BlockBatchStart1 = dcr1->BlockShift(); 
							BlockBatchSize1 = dcr1->size();
							BlockBatchStart2 = 0; 
							dcr2->DumpCompatibleBlocks(cit2->first , Col2 , Row2, dat2, dcr1->MyBlocking, BlockBatchSize2, print_lvl);
							if (BlockBatchSize2 == 0)
							{
								if (!dcr1->MyBlocking.NextBlock())
									break; 
								else 
									continue; 
							}
							if (!dcr1->MyBlocking.back()->NotOver())
								break; 
							
							if (cit1->first.size() == 0)
							{
								dcr1->PartitionAndDumpBlock(cit1->first , Col1 , Row1, dat1);
								Grassman(Row1, dat1, Row2, dat2, BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
											  Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap);								
							}
							else
							{
								dcr1->PartitionAndDumpBlock(cit1->first , Col1 , Row1, dat1);	
								SparseContract(cit1->first, cit2->first, Row1, Col1, dat1, Row2, Col2, dat2,
											   BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
											   Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap);								
							}							
							if ( !dcr1->MyBlocking.NextBlock() )
								break; 	
						}
					}
				}
			}
			else 
			{
				if (print_lvl)
					cout << "Len1 < Len2, picking 2nd branch" << endl;   
				// ReBlock 2 less often because it's larger. 
				for (cit2 = cmap2.begin() ; cit2 != cmap2.end() ; ++cit2)
				{
					dcr2->ReBlock(cit2->first, true);
					if (print_lvl>1)
					{
						cout << "---------Blocking of 2-----------" << endl;  
						dcr2->MyBlocking.print(2);
					}
					
					for (cit1 = cmap1.begin() ; cit1 != cmap1.end() ; ++cit1)
					{
						dcr1->ReBlock(cit1->first, false);
						if (print_lvl>1)
						{
							cout << "---------Blocking of 1-----------" << endl;  
							dcr1->MyBlocking.print(2);
						}
						if (!one.MySym.empty && !two.MySym.empty)
						{
							if (cit2 != cmap2.begin())
							{
								if (one.IsAmp && two.IsAmp)
								{
									if (!one.MySym.CdimsAreCompatible(cit1->first, cit2->first, one.Rank(), two.Rank()))
										continue; 
								}		
								else if (!one.MySym.CdimsAreCompatible(cit1->first, cit2->first, one.MySym.dset, two.MySym.dset))
									continue; 
							}
						}
						//Construct sign,eimap for this pair of contractions
						//Cheap way to pass Dim1 and Dim2
						eff_eimap.clear(); 
						eff_eimap[0] = std::pair<int,int>(one.Rank() , two.Rank() );
						three.MySym.GetSignEimap(eimap , *cit1 , *cit2 , eff_sign , eff_eimap , print_lvl);
						if (print_lvl)
						{
							cout << "CDims: " << eff_sign << endl;  
							printvec(cit1->first);
							printvec(cit2->first);
						} 						
						
						dcr2->MyBlocking.Begin();
						while (true)
						{
							BlockBatchStart2 = dcr2->BlockShift(); 
							BlockBatchSize2 = dcr2->size();
							BlockBatchStart1 = 0; 				
							dcr1->DumpCompatibleBlocks(cit1->first , Col1 , Row1, dat1, dcr2->MyBlocking, BlockBatchSize1, print_lvl);		
							if ( BlockBatchSize1 == 0) 
							{
								if (!dcr2->MyBlocking.NextBlock())
									break; 
								else 
									continue; 
							}
							if (!dcr2->MyBlocking.back()->NotOver())
								break; 							
							if (cit1->first.size() == 0)
							{
								dcr2->PartitionAndDumpBlock(cit2->first , Col2 , Row2, dat2);
								Grassman(Row1, dat1, Row2, dat2, BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
											  Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap);								
							}
							else 
							{
								dcr2->PartitionAndDumpBlock(cit2->first , Col2 , Row2, dat2);	
								SparseContract(cit1->first, cit2->first, Row1, Col1, dat1, Row2, Col2, dat2,
											   BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
											   Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap); 
							}						
							if ( !dcr2->MyBlocking.NextBlock() )
								break; 
						}
					}
				}
			}
			if (print_lvl)
			{
				cout << " ------------------- ContractLayered Result: " << endl; 
				three.print(100);
			}
			// Some of these may 
			
			if (Row1 != NULL) 
				delete [] Row1; 
			if (Col1 != NULL) 
				delete [] Col1; 
			if (Row2 != NULL) 
				delete [] Row2;
			if (Col2 != NULL) 
				delete [] Col2; 
			if (dat1 != NULL) 
				delete [] dat1; 
			if (dat2 != NULL) 
				delete [] dat2;
			if (dcr1 != NULL) 
				delete dcr1; 
			if (dcr2 != NULL) 
				delete dcr2; 
			Row1=NULL; Row2=NULL; Col1=NULL; Col2=NULL; 
			dat1=NULL; dat2=NULL; dcr1=NULL; dcr2=NULL; 
			return; 
		}
		
		// Result is an intermediate
		void ContractLayeredInter(GenTensor& one, GenTensor& two, GenTensor& three , int print_lvl = 0, std::vector<int> DBGDEX = std::vector<int>())
		{
			if (!one.l2norm() || !two.l2norm())
				return;
			// Result must know it's symmetries which may be absent if this was populated then deleted 
			// and re-instantiated. This should be made consistent. 
			three.o_type = o_type;
			three.MySym.AssignPG(osym);
			
			std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >::iterator cit1;
			std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >::iterator cit2;
			map<int,pair<int,int> > eff_eimap; //Effective Eimap. 
			int eff_sign=1;  //Effective sign
			
			DenseLayeredRep* dcr1 = one.GetDLR();
			DenseLayeredRep* dcr2 = two.GetDLR(one.MyTypes);			
			
			if(dcr1 == dcr2)
			{
				cout << "ERROR ONE==TWO in ContractLayeredInter" << endl; 
				system("sleep 100");
			}
			
			double* dat1(NULL);
			double* dat2(NULL);
			int* Col1(NULL);
			int* Row1(NULL);
			int* Col2(NULL);
			int* Row2(NULL);
			dcr1->AllocForPartition((cmap1.begin())->first , Col1, Row1, dat1); 
			dcr2->AllocForPartition((cmap2.begin())->first , Col2, Row2, dat2); 
			
			//Tries to reduce function evalutation. 
			int Rank1 = one.Rank(); 
			int Rank2 = two.Rank(); 
			int Width1 = cmap1.begin()->first.size();
			int Width2 = cmap2.begin()->first.size(); 
			long int ColShift1 = 0; 
			long int ColShift2 = 0; 
			long int RowShift1 = 0; 
			long int RowShift2 = 0; 
			long int BlockBatchStart1 = 0; 
			long int BlockBatchStart2 = 0;
			long int BlockBatchSize1 = 0; 
			long int BlockBatchSize2 = 0;  
			// A bit convoluted, but these temporaries hold the dimensions
			// of One & Two which shouldn't be used for blocking
			std::vector<int> DimUnblockedOne; 
			std::vector<int> DimUnblockedTwo; 
			
			if (dcr1->clength*cmap2.size() > dcr2->clength*cmap1.size())
			{
				//Loop over possible contractions owing to permutational symmetry. 
				for (cit1 = cmap1.begin() ; cit1 != cmap1.end() ; ++cit1)
				{			
					for (cit2 = cmap2.begin() ; cit2 != cmap2.end() ; ++cit2)
					{
						
						bool CdimsCompatible = true; 
						if (!one.MySym.empty && !two.MySym.empty)
						{
							if (cit2 != cmap2.begin())
							{
								if (!one.MySym.OneGroup && !two.MySym.OneGroup)
								{
									if (!one.MySym.CdimsAreCompatible(cit1->first, cit2->first, one.Rank(), two.Rank()))
										CdimsCompatible = false; 
								}		
								else if (!one.MySym.CdimsAreCompatible(cit1->first, cit2->first, one.MySym.dset, two.MySym.dset))
									CdimsCompatible = false; 
							}
						}
						if (!CdimsCompatible)
							continue; 
						eff_eimap.clear(); 
						eff_eimap[0] = std::pair<int,int>(one.Rank() , two.Rank() );
						three.MySym.GetSignEimap(eimap , *cit1 , *cit2 , eff_sign , eff_eimap , 0);
						
						if (print_lvl)
						{
							cout << "CDims: " << eff_sign << endl;  
							printvec(cit1->first);
							printvec(cit2->first);
							printeimap(eff_eimap);
						} 			
						
						if (cit2 == cmap2.begin())
						{
							DimUnblockedOne.clear(); 
							DimUnblockedOne.push_back(one.Rank()); 
							three.MySym.ExternalToArg(eff_eimap, external_eindices, DimUnblockedOne, 1); 
							// Above routine maps to the complement. Gotta map back. 
							one.ComplementaryDims(DimUnblockedOne); 
							dcr1->ReBlock(cit1->first, DimUnblockedOne, true);
						}
						DimUnblockedTwo.clear(); 
						DimUnblockedTwo.push_back(two.Rank());
						three.MySym.ExternalToArg(eff_eimap, external_eindices, DimUnblockedTwo, 2); 
						two.ComplementaryDims(DimUnblockedTwo); 
						dcr2->ReBlock(cit2->first, DimUnblockedTwo, false);
						dcr1->MyBlocking.Begin();
						while (true)
						{
							BlockBatchStart1 = dcr1->BlockShift(); 
							BlockBatchSize1 = dcr1->size();
							BlockBatchStart2 = 0; 
							dcr2->DumpCompatibleBlocks(cit2->first , Col2 , Row2, dat2, dcr1->MyBlocking, BlockBatchSize2 );
							if (BlockBatchSize2 == 0)
							{
								if (!dcr1->MyBlocking.NextBlock())
									break; 
								else 
									continue; 
							}
							if (!dcr1->MyBlocking.back()->NotOver())
								break; 
							if (cit1->first.size() == 0)
							{
								dcr1->PartitionAndDumpBlock(cit1->first , Col1 , Row1, dat1);
								Grassman(Row1, dat1, Row2, dat2, BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
										 Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap, three.MySym, FilterOfOutput);								
							}
							else 
							{
								dcr1->PartitionAndDumpBlock(cit1->first , Col1 , Row1, dat1);	
								SparseContract(cit1->first, cit2->first, Row1, Col1, dat1, Row2, Col2, dat2,
											   BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
											   Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap, 
											   SymmetryOfOutput);
							}
							if (!dcr1->MyBlocking.NextBlock())
								break; 
						}
					}
				}
			}
			else 
			{
				//Loop over possible contractions owing to permutational symmetry. 
				for (cit2 = cmap2.begin() ; cit2 != cmap2.end() ; ++cit2)
				{
					for (cit1 = cmap1.begin() ; cit1 != cmap1.end() ; ++cit1)
					{			
						if (!one.MySym.empty && !two.MySym.empty)
						{
							if (cit1 != cmap1.begin())
							{
								if (one.IsAmp && two.IsAmp)
								{
									if (!one.MySym.CdimsAreCompatible(cit1->first, cit2->first, one.Rank(), two.Rank()))
										continue; 
								}		
								else if (!one.MySym.CdimsAreCompatible(cit1->first, cit2->first, one.MySym.dset, two.MySym.dset))
									continue; 
							}
						}
						
						eff_eimap.clear(); 
						eff_eimap[0] = std::pair<int,int>(one.Rank() , two.Rank() );
						three.MySym.GetSignEimap(eimap , *cit1 , *cit2 , eff_sign , eff_eimap , 0);
						
						if (cit1 == cmap1.begin())
						{
							DimUnblockedTwo.clear(); 
							DimUnblockedTwo.push_back(two.Rank()); 
							three.MySym.ExternalToArg(eff_eimap, external_eindices, DimUnblockedTwo, 2); 
							// Above routine maps to the complement. Gotta map back. 
							two.ComplementaryDims(DimUnblockedTwo); 
							dcr2->ReBlock(cit2->first, DimUnblockedTwo, true);
							if (print_lvl>1)
							{
								cout << "---------Blocking of 2-----------" << endl;  
								dcr2->MyBlocking.print(1);	
							}							
						}
						DimUnblockedOne.clear(); 
						DimUnblockedOne.push_back(one.Rank());
						three.MySym.ExternalToArg(eff_eimap, external_eindices, DimUnblockedOne, 1); 
						one.ComplementaryDims(DimUnblockedOne); 
						dcr1->ReBlock(cit1->first, DimUnblockedOne, false);
						
						dcr2->MyBlocking.Begin();
						while (true)
						{
							BlockBatchStart2 = dcr2->BlockShift(); 
							BlockBatchSize2 = dcr2->size();
							BlockBatchStart1 = 0; 
							dcr1->DumpCompatibleBlocks(cit1->first , Col1 , Row1, dat1, dcr2->MyBlocking, BlockBatchSize1 );
							if (BlockBatchSize1 == 0)
							{
								if (!dcr2->MyBlocking.NextBlock())
									break; 
								else 
									continue; 
							}
							if (!dcr2->MyBlocking.back()->NotOver())
								break; 
							if (cit2->first.size() == 0)
							{
								dcr2->PartitionAndDumpBlock(cit1->first , Col2 , Row2, dat2);
								Grassman(Row1, dat1, Row2, dat2, BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
										 Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap, three.MySym, FilterOfOutput);								
							}
							else 
							{
								dcr2->PartitionAndDumpBlock(cit2->first , Col2 , Row2, dat2);	
								SparseContract(cit1->first, cit2->first, Row1, Col1, dat1, Row2, Col2, dat2,
											   BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
											   Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap, 
											   SymmetryOfOutput);
							}
							if (!dcr2->MyBlocking.NextBlock())
								break; 
						}
					}
				}
			}
			if (print_lvl)
			{
				cout << " ------------------- ContractLayeredInter Result: " << endl; 
				three.print(100);
			}
			if (Row1 != NULL) 
				delete [] Row1; 
			if (Col1 != NULL) 
				delete [] Col1; 
			if (Row2 != NULL) 
				delete [] Row2;
			if (Col2 != NULL) 
				delete [] Col2; 
			if (dat1 != NULL) 
				delete [] dat1; 
			if (dat2 != NULL) 
				delete [] dat2;
			if (dcr1 != NULL) 
				delete dcr1; 
			if (dcr2 != NULL) 
				delete dcr2; 
			Row1=NULL; Row2=NULL; Col1=NULL; Col2=NULL; 
			dat1=NULL; dat2=NULL; dcr1=NULL; dcr2=NULL; 
			return; 
		}

		// Differs from others because it 
		// explicitly Checks for symmetry of Density Matrix. 
		inline void ContractLayeredDM(GenTensor& one, GenTensor& two, GenTensor& three , int print_lvl = 0 , std::vector<int> DBGDEX = std::vector<int>())
		{
			if (!one.l2norm() || !two.l2norm())
				return;
			std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >::iterator cit1;
			std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >::iterator cit2;
			map<int,pair<int,int> > eff_eimap; //Effective Eimap. 
			int eff_sign=1;  //Effective sign
						
			DenseLayeredRep* dcr1 = one.GetDLR();
			DenseLayeredRep* dcr2 = two.GetDLR(one.MyTypes);
			
			if(dcr1 == dcr2)
			{
				cout << "ERROR ONE==TWO in ContractLayeredDM" << endl; 
				system("sleep 100");
			}
			
			double* dat1(NULL); 
			double* dat2(NULL); 
			int* Col1(NULL);
			int* Row1(NULL);
			int* Col2(NULL);
			int* Row2(NULL);
			dcr1->AllocForPartition((cmap1.begin())->first , Col1, Row1, dat1); 
			dcr2->AllocForPartition((cmap2.begin())->first , Col2, Row2, dat2); 
			
			//Tries to reduce function evalutation. 
			int Rank1 = one.Rank(); 
			int Rank2 = two.Rank(); 
			int Width1 = cmap1.begin()->first.size();
			int Width2 = cmap2.begin()->first.size(); 
			long int ColShift1 = 0; 
			long int ColShift2 = 0; 
			long int RowShift1 = 0; 
			long int RowShift2 = 0; 
			long int BlockBatchStart1 = 0; 
			long int BlockBatchStart2 = 0;
			long int BlockBatchSize1 = 0; 
			long int BlockBatchSize2 = 0; 

			if (print_lvl)
			{
				cout << "ContractLayeredDM!" << endl; 
				cout << "Cmap Sizes: " << cmap1.size() << " " << cmap2.size() << endl;  
			}

			//Loop over possible contractions owing to permutational symmetry. 
			for (cit1 = cmap1.begin() ; cit1 != cmap1.end() ; ++cit1)
			{
				dcr1->ReBlock(cit1->first, true);
				for (cit2 = cmap2.begin() ; cit2 != cmap2.end() ; ++cit2)
				{
					if (print_lvl>1)
					{
						cout << "CDims: " << endl;  
						printvec(cit1->first);
						printvec(cit2->first);
					}
					if (!one.MySym.empty && !two.MySym.empty)
						if (cit2 != cmap2.begin())
							if (!one.MySym.CdimsAreCompatible(cit1->first, cit2->first, one.MySym.dset, two.MySym.dset))
								continue; 			
					dcr2->ReBlock(cit2->first, false, print_lvl); // Sorted in DumpCompatibleBlocks
					//Construct sign,eimap for this pair of contractions
					//Cheap way to pass Dim1 and Dim2
					eff_eimap.clear(); 
					eff_eimap[0] = std::pair<int,int>(one.Rank() , two.Rank() );
					three.MySym.GetSignEimap(eimap , *cit1 , *cit2 , eff_sign , eff_eimap , 0);
					dcr1->MyBlocking.Begin();
					while (true)
					{
						BlockBatchStart1 = dcr1->BlockShift(); 
						BlockBatchSize1 = dcr1->size();
						BlockBatchStart2 = 0; 	
						// Where dcr2 is sorted. 
						dcr2->DumpCompatibleBlocks(cit2->first , Col2 , Row2, dat2, dcr1->MyBlocking, BlockBatchSize2, print_lvl);						
						if (BlockBatchSize2 == 0)
						{
							if (!dcr1->MyBlocking.NextBlock())
								break; 
							else 
								continue; 
						}
						if (!dcr1->MyBlocking.back()->NotOver())
							break;
						
						if (cit1->first.size() == 0)
						{
							dcr1->PartitionAndDumpBlock(cit1->first , Col1 , Row1, dat1);
							Grassman(Row1, dat1, Row2, dat2, BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
										  Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap, three.MySym, FilterOfOutput);								
						}
						else 
						{
							// Must sort with the output symmetry && no O/V sorting. 
							dcr1->PartitionAndDumpBlock(cit1->first , Col1 , Row1, dat1);								
							SparseContract(cit1->first, cit2->first, Row1, Col1, dat1, Row2, Col2, dat2,
							   BlockBatchStart1, BlockBatchSize1, BlockBatchStart2, BlockBatchSize2,
							   Rank1 - Width1, Rank2 - Width2, (eff_sign*factor), three, eff_eimap, 
							   three.MySym);
						} 
						if ( !dcr1->MyBlocking.NextBlock() )
							break; 	
					}
				}
			}
			if (Row1 != NULL) 
				delete [] Row1; 
			if (Col1 != NULL) 
				delete [] Col1; 
			if (Row2 != NULL) 
				delete [] Row2;
			if (Col2 != NULL) 
				delete [] Col2; 
			if (dat1 != NULL) 
				delete [] dat1; 
			if (dat2 != NULL) 
				delete [] dat2;
			if (dcr1 != NULL) 
				delete dcr1; 
			if (dcr2 != NULL) 
				delete dcr2; 
			Row1=NULL; Row2=NULL; Col1=NULL; Col2=NULL; 
			dat1=NULL; dat2=NULL; dcr1=NULL; dcr2=NULL; 
			return; 
		}

		// Inner While() loop of Algorithm 4 in the Manuscript. 
		// As implemented by John Parkhill in 2009 at University of California Berkeley. 
		//
		// Requires objects for Dest with just 2 simple functions (SortDex & Increment) 
		// Also CL_DSCAL from BLAS and the C++ std library. 
		// 
		// Performs: Dest += Factor * {Sum}_(Internal1[i] == Internal2[i]) of Tensor1_{External1}^Internal1 * Tensor2_(External2)^Internal2
		// 
		// InternalDims1 is d_1^k from the text. 
		// External1 and Internal1 are the external and internal parts of Tensor1's coordinates. (same order)
		// Tensor1 and Tensor2 are both ordered so that Internal1 and Internal2 are lexicographically ascending. 
		// Data1 are the values at the above coordinate. 
		//
		// Dest is an object for the output tensor with simple increment function x.incr(Coordinate, double Value) incrementing Coordinate with value
		// It also has a Symmetry object MySym which sorts the index over any symmetrical dimensions and returns the sign of that permutation. 
		//
		// eff_eimap gives the source of each result dimension (beginning at 0) ie: 
		// eff_eimap[0] = (0, 2) means the first dimension of output comes from the first tensor's 3rd dimension
		// eff_eimap[1] = (1, 3) means the second dimension of output comes from the second tensor's 4th dimension etc. 
		//
		// Start1, Start2 are offsets for Tensor1 and Tensor2 to support the beginning of a Block. (0 for no Blocks)
		// Sz1, Sz2 are the size of this block (or the size of the whole tensor without blocking).
		// ExtSz1 and ExtSz2 are just the numbers of External dimensions on each tensor. 
		void SparseContract(const std::vector<int>& InternalDims1, const std::vector<int>& InternalDims2,
							int* External1, int* Internal1, double* Data1, int* External2, int* Internal2, double* Data2,
							const long int Start1,const long int Sz1, const long int Start2, const long int Sz2, 
							const int ExtSz1, const int ExtSz2, const double Factor,
							GenTensor& Dest, std::map<int,std::pair<int,int> >& eff_eimap, 
							int print_lvl = 0) const 
		{
			if (Internal1 == NULL || Internal2 == NULL)
				return;
			if (!Sz1 || !Sz2)
				cout << "Arguement Length error in SparseContract" << endl; 
			std::vector<int>::const_iterator DWrite; 
			std::vector<int>::iterator DWrite2; 
			std::vector<int> SourceEDim1(ExtSz1);
			std::vector<int> SourceEDim2(ExtSz2);
			int IntSz1 = InternalDims1.size(); 
			if (ExtSz1)
			{
				DWrite2 = SourceEDim1.begin(); 
				for (int tmp = 0; tmp < (ExtSz1+IntSz1) ; ++tmp)
				{
					if (std::find(InternalDims1.begin(), InternalDims1.end(), tmp) == InternalDims1.end())
					{	
						*DWrite2 = tmp; 
						++DWrite2;
					}
				}		
			}
			if (ExtSz2)
			{
				DWrite2 = SourceEDim2.begin(); 
				for (int tmp = 0; tmp < (ExtSz2+IntSz1) ; ++tmp)
				{
					if (std::find(InternalDims2.begin(), InternalDims2.end(), tmp) == InternalDims2.end())
					{	
						*DWrite2 = tmp; 
						++DWrite2;
					}
				}		
			}
			std::map<int,std::pair<int,int> >::const_iterator MWrite; 
			int TmpInt = 0 ; 
			int TmpInt2 = 0 ; 
			int TmpInt3 = 0 ; 
			std::vector<int> TmpDex(ExtSz1+ExtSz2); 							
			int* Place1 = Internal1 + Start1*IntSz1;
			int* Place2 = Internal2 + Start2*IntSz1;
			double* Work = new double[Sz2]; 
			memset(Work, 0.0, sizeof(double)*Sz2);
			int* CoordWork = new int[(ExtSz1+ExtSz2)*Sz2];  
			memset(CoordWork, 0, sizeof(int)*(ExtSz1+ExtSz2)*Sz2);
			int WorkSize = 0; 
			int* ChunkStart1;
			int* ChunkStart2; 
			int* ChunkEnd1;
			int* ChunkEnd2; 
			while ((Place2 < Internal2+(Start2+Sz2)*IntSz1) && (Place1 < Internal1+(Start1+Sz1)*IntSz1))
			{
				if (std::equal(Place1, Place1+IntSz1, Place2)) 
				{
					if (print_lvl)
					{
						cout << "Contracting Over: " ;
						printvec(std::vector<int>(Place1, Place1+IntSz1), true);
						cout << endl; 
					}
					ChunkStart1 = Place1; 
					ChunkEnd1 = Place1; 
					ChunkStart2 = Place2; 
					ChunkEnd2 = Place2; 					
					// Find the end of this Chunk... 
					while( (ChunkEnd2 < Internal2+(Start2+Sz2)*IntSz1) && std::equal(Place1, Place1+IntSz1, ChunkEnd2))
						ChunkEnd2 += IntSz1;
					while( (ChunkEnd1 < Internal1+(Start1+Sz1)*IntSz1) && std::equal(Place1, Place1+IntSz1, ChunkEnd1))
						ChunkEnd1 += IntSz1;
					int AbsoluteStart1 = (ChunkStart1 - Internal1)/IntSz1;
					int AbsoluteStart2 = (ChunkStart2 - Internal2)/IntSz1;
					int AbsoluteEnd1 = (ChunkEnd1 - Internal1)/IntSz1 - 1;
					int AbsoluteEnd2 = (ChunkEnd2 - Internal2)/IntSz1 - 1;
					
					if (print_lvl > 2)
					{
						cout << "Starts, End: " << AbsoluteStart1 << " : " << AbsoluteEnd1  <<  " :: " << AbsoluteStart2 << " : " << AbsoluteEnd2 << endl; 
						for (int SpotIn1 = AbsoluteStart1 ; SpotIn1 <= AbsoluteEnd1 ; ++SpotIn1)
						{
							cout << "Int1: ";
							for (int tm = 0 ; tm < IntSz1 ; ++tm)
								cout << " " << Internal1[SpotIn1*IntSz1+tm]; 
							cout << " Ext1: ";
							for (int tm = 0 ; tm < ExtSz1 ; ++tm)
								cout << " " << External1[SpotIn1*ExtSz1+tm];
							cout << endl; 
						}
						for (int SpotIn1 = AbsoluteStart2 ; SpotIn1 <= AbsoluteEnd2 ; ++SpotIn1)
						{
							cout << "Int2: ";
							for (int tm = 0 ; tm < IntSz1 ; ++tm)
								cout << " " << Internal2[SpotIn1*IntSz1+tm]; 
							cout << " Ext2: ";
							for (int tm = 0 ; tm < ExtSz1 ; ++tm)
								cout << " " << External2[SpotIn1*ExtSz2+tm];
							cout << endl; 
						}
					}
					int ChunkSize2 = AbsoluteEnd2 - AbsoluteStart2 + 1;
					for (MWrite = eff_eimap.begin(); MWrite != eff_eimap.end(); ++MWrite)
					{
						if (MWrite->second.first > 0)
						{
							TmpInt = MWrite->first;
							TmpInt2 = (std::find(SourceEDim2.begin(), SourceEDim2.end(), MWrite->second.second) - SourceEDim2.begin());
							for (int Spot2x = 0 ; Spot2x < ChunkSize2 ; ++Spot2x)
								CoordWork[Spot2x*(ExtSz1+ExtSz2)+TmpInt] = External2[(Spot2x+AbsoluteStart2)*ExtSz2+TmpInt2]; 
						}
					}					
					for (int SpotIn1 = AbsoluteStart1 ; SpotIn1 <= AbsoluteEnd1 ; ++SpotIn1)
					{
						if (fabs(Data1[SpotIn1]) < pow(10.0,-17.0))
							continue; 
						//Make 1's Index Contribution. 
						for (MWrite = eff_eimap.begin(); MWrite != eff_eimap.end(); ++MWrite)
						{
							if (MWrite->second.first == 0)
							{
								TmpInt = MWrite->first; 
								TmpInt2 = (std::find(SourceEDim1.begin(), SourceEDim1.end(), MWrite->second.second) - SourceEDim1.begin());
								TmpInt3 = External1[SpotIn1*ExtSz1+TmpInt2];
								for (int Spot2 = 0 ; Spot2 < ChunkSize2 ; ++Spot2)
									CoordWork[Spot2*(ExtSz1+ExtSz2)+TmpInt] = TmpInt3; 
							}
						}
						memcpy(Work, Data2+AbsoluteStart2, ChunkSize2*sizeof(double));
						CL_DSCAL(ChunkSize2, Factor*Data1[SpotIn1], Work, 1);
						for (int Spot2 = 0 ; Spot2 < ChunkSize2 ; ++Spot2)
						{
							if (fabs(Work[Spot2]) < pow(10.0,-15.0))
								continue; 
							TmpDex.assign(CoordWork+Spot2*(ExtSz1+ExtSz2), CoordWork+(Spot2+1)*(ExtSz1+ExtSz2)); 
							Work[Spot2] *= Dest.MySym.SortDex(TmpDex); 
							if (fabs(Work[Spot2]) > pow(10.0,-15.0))
							{
								if (FilterOfOutput.obeys(TmpDex))
									Dest.incr(TmpDex, Work[Spot2]);
							}
						}
					}
					// Place1 & 2 are moved past the end of this chunk
					Place1 = ChunkEnd1;
					Place2 = ChunkEnd2;
				}
				else 
				{
					if (std::lexicographical_compare(Place1, Place1+IntSz1, Place2,Place2+IntSz1))
						Place1+=IntSz1;
					else
						Place2+=IntSz1; 
				}
			}
			delete [] CoordWork; 
			delete [] Work; 
		}
		
		// Sparse Contraction of two vectorized tensors
		// blocking shift is present. 
		void SparseContract(const std::vector<int>& InternalDims1, const std::vector<int>& InternalDims2,
							int* External1, int* Internal1, double* Data1, int* External2, int* Internal2, double* Data2,
							const long int Start1,const long int Sz1, const long int Start2, const long int Sz2, 
							const int ExtSz1, const int ExtSz2, const double Factor, GenTensor& Dest,
						    std::map<int,std::pair<int,int> >& eff_eimap, TensorSymmetry& OutSym, 
							int print_lvl = 0) const 
		{
			if (Internal1 == NULL || Internal2 == NULL)
				return;
			if (!Sz1 || !Sz2)
				cout << "Arguement Length error in SparseContract" << endl; 
			std::vector<int>::const_iterator DWrite; 
			std::vector<int>::iterator DWrite2; 
			std::vector<int> SourceEDim1(ExtSz1);
			std::vector<int> SourceEDim2(ExtSz2);
			int IntSz1 = InternalDims1.size(); 
			if (ExtSz1)
			{
				DWrite2 = SourceEDim1.begin(); 
				for (int tmp = 0; tmp < (ExtSz1+IntSz1) ; ++tmp)
				{
					if (std::find(InternalDims1.begin(), InternalDims1.end(), tmp) == InternalDims1.end())
					{	
						*DWrite2 = tmp; 
						++DWrite2;
					}
				}		
			}
			if (ExtSz2)
			{
				DWrite2 = SourceEDim2.begin(); 
				for (int tmp = 0; tmp < (ExtSz2+IntSz1) ; ++tmp)
				{
					if (std::find(InternalDims2.begin(), InternalDims2.end(), tmp) == InternalDims2.end())
					{	
						*DWrite2 = tmp; 
						++DWrite2;
					}
				}		
			}
			std::map<int,std::pair<int,int> >::const_iterator MWrite; 
			int TmpInt = 0 ; 
			int TmpInt2 = 0 ; 
			int TmpInt3 = 0 ; 
			std::vector<int> TmpDex(ExtSz1+ExtSz2); 							
			int* Place1 = Internal1 + Start1*IntSz1;
			int* Place2 = Internal2 + Start2*IntSz1;
			double* Work = new double[Sz2]; 
			memset(Work, 0.0, sizeof(double)*Sz2);
			int* CoordWork = new int[(ExtSz1+ExtSz2)*Sz2];  
			memset(CoordWork, 0, sizeof(int)*(ExtSz1+ExtSz2)*Sz2);
			int WorkSize = 0; 
			int* ChunkStart1;
			int* ChunkStart2; 
			int* ChunkEnd1;
			int* ChunkEnd2; 
			while ((Place2 < Internal2+(Start2+Sz2)*IntSz1) && (Place1 < Internal1+(Start1+Sz1)*IntSz1))
			{
				if (std::equal(Place1, Place1+IntSz1, Place2)) 
				{
					if (print_lvl>1)
					{
						cout << "Contracting Over: " ;
						printvec(std::vector<int>(Place1, Place1+IntSz1), true);
						cout << endl; 
					}
					ChunkStart1 = Place1; 
					ChunkEnd1 = Place1; 
					ChunkStart2 = Place2; 
					ChunkEnd2 = Place2; 					
					// Find the end of this Chunk... 
					while( (ChunkEnd2 < Internal2+(Start2+Sz2)*IntSz1) && std::equal(Place1, Place1+IntSz1, ChunkEnd2))
						ChunkEnd2 += IntSz1;
					while( (ChunkEnd1 < Internal1+(Start1+Sz1)*IntSz1) && std::equal(Place1, Place1+IntSz1, ChunkEnd1))
						ChunkEnd1 += IntSz1;
					int AbsoluteStart1 = (ChunkStart1 - Internal1)/IntSz1;
					int AbsoluteStart2 = (ChunkStart2 - Internal2)/IntSz1;
					int AbsoluteEnd1 = (ChunkEnd1 - Internal1)/IntSz1 - 1;
					int AbsoluteEnd2 = (ChunkEnd2 - Internal2)/IntSz1 - 1;
					
					if (print_lvl>2)
					{
						cout << "Starts, End: " << AbsoluteStart1 << " : " << AbsoluteEnd1  <<  " :: " << AbsoluteStart2 << " : " << AbsoluteEnd2 << endl; 
						for (int SpotIn1 = AbsoluteStart1 ; SpotIn1 <= AbsoluteEnd1 ; ++SpotIn1)
						{
							cout << "Int1: ";
							for (int tm = 0 ; tm < IntSz1 ; ++tm)
								cout << " " << Internal1[SpotIn1*IntSz1+tm]; 
							cout << " Ext1: ";
							for (int tm = 0 ; tm < ExtSz1 ; ++tm)
								cout << " " << External1[SpotIn1*ExtSz1+tm];
							cout << endl; 
						}
						for (int SpotIn1 = AbsoluteStart2 ; SpotIn1 <= AbsoluteEnd2 ; ++SpotIn1)
						{
							cout << "Int2: ";
							for (int tm = 0 ; tm < IntSz1 ; ++tm)
								cout << " " << Internal2[SpotIn1*IntSz1+tm]; 
							cout << " Ext2: ";
							for (int tm = 0 ; tm < ExtSz1 ; ++tm)
								cout << " " << External2[SpotIn1*ExtSz2+tm];
							cout << endl; 
						}
					}
					int ChunkSize2 = AbsoluteEnd2 - AbsoluteStart2 + 1;
					for (MWrite = eff_eimap.begin(); MWrite != eff_eimap.end(); ++MWrite)
					{
						if (MWrite->second.first > 0)
						{
							TmpInt = MWrite->first;
							TmpInt2 = (std::find(SourceEDim2.begin(), SourceEDim2.end(), MWrite->second.second) - SourceEDim2.begin());
							for (int Spot2x = 0 ; Spot2x < ChunkSize2 ; ++Spot2x)
								CoordWork[Spot2x*(ExtSz1+ExtSz2)+TmpInt] = External2[(Spot2x+AbsoluteStart2)*ExtSz2+TmpInt2]; 
						}
					}					
					for (int SpotIn1 = AbsoluteStart1 ; SpotIn1 <= AbsoluteEnd1 ; ++SpotIn1)
					{
						if (fabs(Data1[SpotIn1]) < pow(10.0,-17.0))
							continue; 
						//Make 1's Index Contribution. 
						for (MWrite = eff_eimap.begin(); MWrite != eff_eimap.end(); ++MWrite)
						{
							if (MWrite->second.first == 0)
							{
								TmpInt = MWrite->first; 
								TmpInt2 = (std::find(SourceEDim1.begin(), SourceEDim1.end(), MWrite->second.second) - SourceEDim1.begin());
								TmpInt3 = External1[SpotIn1*ExtSz1+TmpInt2];
								for (int Spot2 = 0 ; Spot2 < ChunkSize2 ; ++Spot2)
									CoordWork[Spot2*(ExtSz1+ExtSz2)+TmpInt] = TmpInt3; 
							}
						}
						memcpy(Work, Data2+AbsoluteStart2, ChunkSize2*sizeof(double));
						CL_DSCAL(ChunkSize2, Factor*Data1[SpotIn1], Work, 1);
						for (int Spot2 = 0 ; Spot2 < ChunkSize2 ; ++Spot2)
						{
							if (fabs(Work[Spot2]) < pow(10.0,-15.0))
								continue; 
							TmpDex.assign(CoordWork+Spot2*(ExtSz1+ExtSz2), CoordWork+(Spot2+1)*(ExtSz1+ExtSz2)); 
							Work[Spot2] *= OutSym.SortDex(TmpDex, OutSym.dset);
							if (fabs(Work[Spot2]) > pow(10.0,-15.0))
							{
								if (FilterOfOutput.obeys(TmpDex))
									Dest.incr(TmpDex, Work[Spot2]);
							}
						}
					}
					// Place1 & 2 are moved past the end of this chunk
					Place1 = ChunkEnd1;
					Place2 = ChunkEnd2;
				}
				else 
				{
					if (std::lexicographical_compare(Place1, Place1+IntSz1, Place2,Place2+IntSz1))
						Place1+=IntSz1;
					else
						Place2+=IntSz1; 
				}
			}
			delete [] CoordWork; 
			delete [] Work; 
		}		
		
		// Does the work of a grassman product
		// (Contraction with no common indices)
		void Grassman( int* External1, double* Data1, int* External2, double* Data2,
			const long int Start1,const long int Sz1, const long int Start2,const long int Sz2, 
			const int ExtSz1, const int ExtSz2, const double Factor, GenTensor& Dest, 
			std::map<int,std::pair<int,int> >& eff_eimap)
		{		
			std::vector<int>::const_iterator DWrite; 
			std::vector<int>::iterator DWrite2; 
			std::vector<int> SourceEDim1(ExtSz1);
			std::vector<int> SourceEDim2(ExtSz2);
			int IntSz1 = 0; 
			if (ExtSz1)
			{
				DWrite2 = SourceEDim1.begin(); 
				for (int tmp = 0; tmp <(ExtSz1+IntSz1) ; ++tmp)
				{
					if (true)
					{	
						*DWrite2 = tmp; 
						++DWrite2;
					}
				}		
			}
			if (ExtSz2)
			{
				DWrite2 = SourceEDim2.begin(); 
				for (int tmp = 0; tmp <(ExtSz2+IntSz1) ; ++tmp)
				{
					if (true)
					{	
						*DWrite2 = tmp; 
						++DWrite2;
					}
				}		
			}
			std::map<int,std::pair<int,int> >::const_iterator MWrite; 
			int TmpInt = 0 ; 
			int TmpInt2 = 0 ; 
			int TmpInt3 = 0 ; 
			std::vector<int> TmpDex(ExtSz1+ExtSz2); 
			double* Work = new double[Sz2]; 
			int* CoordWork = new int[(ExtSz1+ExtSz2)*Sz2];  
			for (int Spot1 = 0 ; Spot1 < Sz1 ; ++Spot1)
			{
				if (Data1[Start1+Spot1] == 0)
					continue; 
				//Make 1's Index Contribution. 
				for (MWrite = eff_eimap.begin(); MWrite != eff_eimap.end(); ++MWrite)
				{
					if (MWrite->second.first == 0)
					{
						TmpInt = MWrite->first; 
						TmpInt2 = (std::find(SourceEDim1.begin(), SourceEDim1.end(), MWrite->second.second) - SourceEDim1.begin());
						TmpInt3 = External1[(Spot1+Start1)*ExtSz1+TmpInt2];
						for (int Spot2 = 0 ; Spot2 < Sz2 ; ++Spot2)
							CoordWork[Spot2*(ExtSz1+ExtSz2)+TmpInt] = TmpInt3; 
					}
					else 
					{
						TmpInt = MWrite->first;
						TmpInt2 = (std::find(SourceEDim2.begin(), SourceEDim2.end(), MWrite->second.second) - SourceEDim2.begin());
						for (int Spot2x = 0 ; Spot2x < Sz2 ; ++Spot2x)
						{
							TmpInt3 = External2[(Spot2x+Start2)*ExtSz2+TmpInt2];
							CoordWork[Spot2x*(ExtSz1+ExtSz2)+TmpInt] = TmpInt3; 
						}
					}
				}
				memcpy(Work, Data2+Start2, Sz2*sizeof(double));
				CL_DSCAL(Sz2, Factor*Data1[Start1+Spot1], Work, 1);
				for (int Spot2 = 0 ; Spot2 < Sz2 ; ++Spot2)
				{
					if(Work[Spot2] != 0)
					{
						TmpDex.assign(CoordWork+Spot2*(ExtSz1+ExtSz2), CoordWork+(Spot2+1)*(ExtSz1+ExtSz2)); 
						Work[Spot2] *= Dest.MySym.SortDex(TmpDex); 
						Dest.incr(TmpDex, Work[Spot2]);
					}
				} 
			}
			delete [] CoordWork; 
			delete [] Work; 
		}

		// Does the work of a grassman product
		// (Contraction with no common indices)
		void Grassman( int* External1, double* Data1, int* External2, double* Data2,
			const long int Start1,const long int Sz1, const long int Start2,const long int Sz2, 
			const int ExtSz1, const int ExtSz2, const double Factor, GenTensor& Dest, 
			std::map<int,std::pair<int,int> >& eff_eimap, TensorSymmetry& OutSym, const Filters& OutFilter)
		{
			std::vector<int>::const_iterator DWrite; 
			std::vector<int>::iterator DWrite2; 
			std::vector<int> SourceEDim1(ExtSz1);
			std::vector<int> SourceEDim2(ExtSz2);
			int IntSz1 = 0; 
			if (ExtSz1)
			{
				DWrite2 = SourceEDim1.begin(); 
				for (int tmp = 0; tmp <(ExtSz1+IntSz1) ; ++tmp)
				{
					if (true)
					{	
						*DWrite2 = tmp; 
						++DWrite2;
					}
				}		
			}
			if (ExtSz2)
			{
				DWrite2 = SourceEDim2.begin(); 
				for (int tmp = 0; tmp <(ExtSz2+IntSz1) ; ++tmp)
				{
					if (true)
					{	
						*DWrite2 = tmp; 
						++DWrite2;
					}
				}		
			}
			std::map<int,std::pair<int,int> >::const_iterator MWrite; 
			int TmpInt = 0 ; 
			int TmpInt2 = 0 ; 
			int TmpInt3 = 0 ; 
			std::vector<int> TmpDex(ExtSz1+ExtSz2); 
			double* Work = new double[Sz2]; 
			int* CoordWork = new int[(ExtSz1+ExtSz2)*Sz2];  
			for (int Spot1 = 0 ; Spot1 < Sz1 ; ++Spot1)
			{
				if (Data1[Start1+Spot1] == 0)
					continue; 
				//Make 1's Index Contribution. 
				for (MWrite = eff_eimap.begin(); MWrite != eff_eimap.end(); ++MWrite)
				{
					if (MWrite->second.first == 0)
					{
						TmpInt = MWrite->first; 
						TmpInt2 = (std::find(SourceEDim1.begin(), SourceEDim1.end(), MWrite->second.second) - SourceEDim1.begin());
						TmpInt3 = External1[(Spot1+Start1)*ExtSz1+TmpInt2];
						for (int Spot2 = 0 ; Spot2 < Sz2 ; ++Spot2)
							CoordWork[Spot2*(ExtSz1+ExtSz2)+TmpInt] = TmpInt3; 
					}
					else 
					{
						TmpInt = MWrite->first;
						TmpInt2 = (std::find(SourceEDim2.begin(), SourceEDim2.end(), MWrite->second.second) - SourceEDim2.begin());
						for (int Spot2x = 0 ; Spot2x < Sz2 ; ++Spot2x)
						{
							TmpInt3 = External2[(Spot2x+Start2)*ExtSz2+TmpInt2];
							CoordWork[Spot2x*(ExtSz1+ExtSz2)+TmpInt] = TmpInt3; 
						}
					}
				}
				memcpy(Work, Data2+Start2, Sz2*sizeof(double));
				CL_DSCAL(Sz2, Factor*Data1[Start1+Spot1], Work, 1);
				for (int Spot2 = 0 ; Spot2 < Sz2 ; ++Spot2)
				{
					if(Work[Spot2] != 0)
					{
						TmpDex.assign(CoordWork+Spot2*(ExtSz1+ExtSz2), CoordWork+(Spot2+1)*(ExtSz1+ExtSz2)); 
						Work[Spot2] *= OutSym.SortDex(TmpDex, OutSym.dset); 
						if(Work[Spot2] != 0)
						{
							if (OutFilter.obeys(TmpDex))
								Dest.incr(TmpDex, Work[Spot2]);
						}
					}
				} 
			}
			delete [] CoordWork; 
			delete [] Work; 
		}
		
		// Very specialized contraction algorithm for problems of the sort: 
		// T(abcdef)T(abcd) += R(de)  {summation implied}
		// if arguements have symmetry these symmetries MUST BE THE SAME. 
		// result must be without symmetry, orbital-type. etc. 
		void DotContract(GenTensor& one, GenTensor& two, GenTensor& three, double fact = 1.0) 
		{
			if (one.l2norm() == 0.0 || two.l2norm() == 0.0)
				return; 
			if (one.Rank() < two.Rank())
				cout << "Rank mismatch 1 in DotContract... " << endl; 
			if (one.Rank() - two.Rank() != three.Rank())
				cout << "Rank mismatch 2 in DotContract... " << endl; 
			
//			cout << "DotContract: " << one.l2norm() << " " << two.l2norm() << endl;  
//			cout << "Factor: " << one.MySym.order()*fact << endl;  			

			int srcrnk = two.Rank(); 
			std::vector<int> i1(one.Rank()); std::vector<int> i2(two.Rank()); std::vector<int> i3(three.Rank()); 
			double pgfac = one.MySym.order(); 
			one.begin(); 
			while (true)
			{
				i1 = one.ind(); 
				std::copy(i1.begin(), i1.begin()+srcrnk , i2.begin()); 
				if (two.nonzero(i2))
				{
					std::copy(i1.begin()+srcrnk, i1.end() , i3.begin()); 	
				//	printvec(i3,1); cout << " " << one.val()*two.val(i2)*pgfac << endl;   
					three.incr(i3,one.val()*two.val(i2)*pgfac*fact);
				}
				if (!one.I())
					break; 
			}
		};
		
		// Debug routines follow ----------------------------------
		inline void printvec(const std::vector<int> arg, bool IfNoLine = false) const
		{	
			std::vector<int>::const_iterator pI2; 
			for (pI2 = arg.begin() ; pI2 != arg.end() ; ++pI2)
			{
				cout << " " << *pI2; 
			}
			if (!IfNoLine)
				cout << endl; 
		}
		inline void printcmap(const std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >& arg) const
		{
			std::vector<std::pair<std::vector<int>,std::pair<std::vector<int> ,std::vector<int> > > >::const_iterator MWrite; 
			cout << "Cmap: " << endl;  
			for (MWrite = arg.begin(); MWrite != arg.end(); ++MWrite)
			{
				printvec(MWrite->first,true); 
				cout << " : ";
				printvec(MWrite->second.first,true); 
				cout << " , " ; 
				printvec(MWrite->second.second,true);
				cout << endl;   
			}
		}
		void printcmaps()
		{
			printcmap(cmap1); 
			printcmap(cmap2); 
		} 
		inline void printeimap(const std::map<int,std::pair<int,int> >& arg) const
		{
			std::map<int,std::pair<int,int> >::const_iterator MWrite; 
			cout << "Eimap: " << endl;  
			for (MWrite = arg.begin(); MWrite != arg.end(); ++MWrite)
			{
				cout << MWrite->first << " : " << MWrite->second.first << " , " << MWrite->second.second << endl;   
			}
		}
		inline void printsov(std::set<std::vector<int> > arg)
		{	
			std::set<std::vector<int> >::const_iterator pI1; 
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
		inline void printsov(std::list<std::vector<int> > arg)
		{	
			std::list<std::vector<int> >::const_iterator pI1; 
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
	};

// Can Store expanded contraction maps and
// restore previously computed data. 
class Contractions : public std::map<int,Contraction>
{
public:
	std::string MyName; 
	Contractions() 
	{}
	Contractions(const char* arg) : MyName(arg) 
	{}
	void Name(char* Nm)
	{
		MyName = std::string(Nm);
	}
	void Write(Streamer& Dest)
	{ 
		Dest.WriteToStream("Contractions: ");
		Dest.WriteToStream((int)size());
		Dest.WriteToStream(" { ");
		for (Contractions::iterator iter = begin(); iter != end() ; ++iter)
		{
			iter->second.Write(Dest);
			Dest.Flush(); 
		}
		Dest.WriteToStream(" } ");
		Dest.Flush();
	}
	void WriteIfAbsent()
	{
		Streamer Tmp; 
		if (Tmp.OpenIfAbsent(MyName.c_str()))
		{
			Write(Tmp); 
			// debug 
/*			cout << "Comparing Written with Current... " << endl; 
			Contractions Hold(*this);
			clear(); 
			ReadIfExists(); 
			Compare(Hold);
*/			//\debug
		}
		else 
			return; 
	}
	void Compare(Contractions& other)
	{
		Contractions::iterator iter; 
		for (iter = other.begin(); iter != other.end();++iter)
		{
			if (!(iter->second.cmap1 == ((*this)[iter->first]).cmap1))
			{
				cout << "Found Difference (1)... " << iter->first << endl; 
				iter->second.cmap1.print();
				cout << " --- " << endl; 
				((*this)[iter->first]).cmap1.print(); 
				system("sleep 10"); 
			}
			if (!(iter->second.cmap2 == ((*this)[iter->first]).cmap2))
			{
				cout << "Found Difference (1)... " << iter->first << endl; 
				iter->second.cmap2.print();
				cout << " --- " << endl; 
				((*this)[iter->first]).cmap2.print(); 
				system("sleep 10"); 
			}
			
		}
	}
	void ReadIfExists()
	{
		if (!empty())
			return; 
		Streamer Tmp; 
		if (Tmp.OpenIfExists(MyName.c_str()))
			Read(Tmp); 
		else 
			return; 
	}
	void Read(Streamer& Source)
	{
		clear(); 
		ifstream* is = Source.is; 
		if (!(is->good()))
		{
			cout << "Bad Source Contractions::Read()" << endl;  
			return; 
		}
		string A = Source.ReadNextString(); 
		string B("Contractions:"); 
		if (A != B)
		{
			cout << "File Doesn't contain Contractions... " << endl; 
			cout << A << endl; 
			cout << B << endl; 
		}
		int Sz = Source.ReadNextInt(); 
		// Read the next {
		A = Source.ReadNextString(); 
		for (int TermNum = 1 ; TermNum <= Sz ;++TermNum)
		{
			// format is {contraction #, cmap1, cmap2} ...
			(*this)[TermNum].Read(Source);
		}
	}
};

template<class mbt> double DivideAndSwap(mbt& oldt, mbt& newt, mbt& denom) 
{
	newt /= denom;
	oldt -= newt;
	double toret = oldt.l2norm(); 
	oldt = newt; 
	return toret; 
};
template<class mbt> double DivideAndSwap(mbt& oldt, mbt& newt, mbt& denom,mbt& error,int diis_type) 
{
	error.set2zero(); 
	newt /= denom;
	oldt -= newt;
	double toret = oldt.l2norm(); 
	if( 2 == diis_type ) //For diis2 need sqrt(dOOVV)
    {
		mbt sqrt_dOV(denom);
		sqrt_dOV.Sqrt();
		oldt *= sqrt_dOV; 
    }
	error = oldt; 
	oldt = newt; 
	return toret; 
};
template<class mbt> void DivideAndSwapBT(mbt& oldt, mbt& newt, mbt& denom,mbt& error,int diis_type) 
{
	error.Set(); 
	newt /= denom;
	oldt -= newt;
	if( 2 == diis_type ) //For diis2 need sqrt(dOOVV)
    {
		mbt sqrt_dOV(denom);
		sqrt_dOV.CopyFrom(denom);
		sqrt_dOV.Sqrt();
		oldt.VectorMultiply(sqrt_dOV); 
    }
	error = oldt; 
	oldt = newt; 
};

////////////////////////////////////////////////////////////////////////////////
// Set the result tensor as the sum of two lower rank tensors ie:
// tijab = (tia) + (tib)  
// First arguement must have greater rank. 

template<class mbt1,class mbt2,class mbt3> void TensorSum(mbt1& one,mbt2& two,mbt3& result, double thresh = 0.0) 
{
	int num_filtered=0; 
	result.set2zero();
	
	if (result.physicalLength() < 1)
	{
		cout << "*** No Structure in TensorSum *** " << endl;  	
		system("sleep 100"); 
	}
	else // If the result structure is already set it's easy to get the tensor product. 
	{
		double filter = 0.0;
		typename std::vector<int> p1(one.Rank(),0);
		typename std::vector<int> p2(two.Rank(),0);
		result.begin(); 
		while(true)
		{
			int placein1 = 0;
			int placein2 = 0;
			for(int iter = 0;iter < result.Rank() ; iter++ )
			{
				if (iter < (mbt3::Dimensions)/2) // Dim is O.
				{
					if (iter%2 == 0)
					{
						p1[placein1] = result.ind()[iter];
						placein1++;
					}	
					else
					{
						p2[placein2] = result.ind()[iter];
						placein2++;
					}
				}
				else // Dim is V.
				{
					int iter2 = (iter - ((result.Rank())/2));
					if (iter2%2 == 0)
					{
						p1[placein1] = result.ind()[iter];
						placein1++;
					}	
					else
					{
						p2[placein2] = result.ind()[iter];
						placein2++;
					}
				}
			}
			filter = one.val(p1) + two.val(p2);
			if (fabs(filter) >= thresh)
				result.set(filter);
			else 
				result.set(thresh);
			if (!result.I())
				break; 
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Set the result tensor as the direct sum of several copies of just one matrix. 
// tijab = (tia) + (tjb) 
// tijkabc = (tia) + (tjb) + (tkc)
// Result tensor is symmetric. 
// Cost is size of result tensor. 

template<class mbt1,class mbt3> void TensorSum(mbt1& one,mbt3& result, double thresh = 0.0) 
{
	int num_filtered=0; 
	result.set2zero();
	if(0 != (result.Rank())%2)
		cout << "Rank mismatch in TensorProduct" << endl;  
	
	//Number of multiples needed
	int multiples = (result.Rank())/2;
	
	if (result.physicalLength() < 1)
	{
		cout << "This Version of TensorSum Requires structure... " << endl; 
	}
	else 
	{
		double filter = 0.0;
		std::map<int,std::vector<int> > images;  //Indices this element will come from		
		
		result.begin(); 
		while(true)
		{
			images.clear(); 
			for(int iter = 0;iter < mbt3::Dimensions ; iter++ )
				if (iter < (result.Rank())/2) // Index is O.
					images[iter%multiples] = std::vector<int>(one.Rank(),0);
			
			for(int iter = 0;iter < mbt3::Dimensions ; iter++ )
			{
				if (iter < (result.Rank())/2) // Index is O.
				{
					(images[iter%multiples])[0] = result.ind()[iter];
				}
				else //Index is V.
				{
					int iter2 = (iter - ((result.Rank())/2));
					(images[iter2%multiples])[1] = result.ind()[iter];				
				}
			}
			filter = 0.0;
			for (std::map<int,std::vector<int> >::iterator iter = images.begin(); iter != images.end() ; iter++)
			{
				filter += one.val(iter->second);
			}
			if (fabs(filter) >= thresh)
				result.set(filter);
			else 
				result.set(thresh);
			if (!result.I())
				break; 
		}
	}
}

#endif
#endif
