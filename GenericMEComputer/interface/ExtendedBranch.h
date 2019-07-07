/** \class ExtendedBranch
*
*
*  \author N. Amapane - Torino
*  \author U. Sarica - JHU
*/
#ifndef EXTENDEDBRANCH_H
#define EXTENDEDBRANCH_H

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cassert>
#include "TString.h"
#include "TTree.h"


namespace BranchHelpers{
  enum BranchTypes{
    bBool,
    bUChar, bChar,
    bUShort, bShort,
    bUInt, bInt,
    // No types are implemented for unsigned long or long because they do not contain a portable implementation in ROOT. 
    // See note in the RtypesCore.h ROOT header file.
    bULongLong, bLongLong,
    bFloat,
    bDouble,
    nBranchTypes
  };

  // Type checking, adaptation of suggestion by Michael Aaron Safyan at http://stackoverflow.com/questions/2728078/is-there-an-easy-way-to-check-a-fundamental-type
  template<BranchTypes VAL> struct constant_type_value{
    static const BranchTypes type = VAL;
  };
  typedef constant_type_value<nBranchTypes> no_type;
  typedef constant_type_value<bBool> bool_type;
  typedef constant_type_value<bUChar> uchar_type;
  typedef constant_type_value<bChar> char_type;
  typedef constant_type_value<bUShort> ushort_type;
  typedef constant_type_value<bShort> short_type;
  typedef constant_type_value<bUInt> uint_type;
  typedef constant_type_value<bInt> int_type;
  typedef constant_type_value<bULongLong> ulonglong_type;
  typedef constant_type_value<bLongLong> longlong_type;
  typedef constant_type_value<bFloat> float_type;
  typedef constant_type_value<bDouble> double_type;

  template<typename T> struct hasType : public no_type{};
  template<> struct hasType<Bool_t> : public bool_type{};
  template<> struct hasType<UChar_t> : public uchar_type{};
  template<> struct hasType<Char_t> : public char_type{};
  template<> struct hasType<UShort_t> : public ushort_type{};
  template<> struct hasType<Short_t> : public short_type{};
  template<> struct hasType<UInt_t> : public uint_type{};
  template<> struct hasType<Int_t> : public int_type{};
  template<> struct hasType<ULong64_t> : public ulonglong_type{};
  template<> struct hasType<Long64_t> : public longlong_type{};
  template<> struct hasType<Float_t> : public float_type{};
  template<> struct hasType<Double_t> : public double_type{};

  template<typename varType> class ExtendedBranch{

  protected:

    TTree* theTree;
    varType* defVal;

  public:

    TString bname;
    varType value;
    std::vector<varType> valueArray;
    BranchHelpers::BranchTypes btype;
    Bool_t isVector;

  public:

    void setValue(varType inVal){
      if (isVector) valueArray.push_back(inVal);
      else value = inVal;
    }
    varType getDefVal() const{
      if (!isVector && defVal) return *defVal;
      else return 0;
    }
    varType getVal() const{
      if (!isVector) return value;
      else if (!valueArray.empty()) return valueArray.back();
      else return 0;
    }
    std::vector<varType> const& getArray() const{ return valueArray; }
    void reset(){
      valueArray.clear();
      value = this->getDefVal();
    }

    TBranch* getBranch(){ return theTree->GetBranch(bname); }

    virtual void Print(){
      std::cout << "**********\n";
      std::cout << "ExtendedBranch summary for " << bname << ":\n";
      std::cout << "\tType=" << btype << "\n";
      if (isVector){
        std::cout << "\tValues = ";
        for (varType const& tmp_val:valueArray) std::cout << tmp_val << " ";
        std::cout << "\n";
      }
      else std::cout << "\tValue = " << value << " (default = " << (defVal ? *defVal : 0) << ")\n";
      std::cout << "\tTree address: " << theTree << "\n";
      std::cout << "**********";
      std::cout << std::endl;
    }

    ExtendedBranch(TTree* theTree_, TString bname_, varType const& defVal_, Bool_t isVector_=false) :
      theTree(theTree_), defVal(nullptr), bname(bname_), isVector(isVector_)
    {
      assignBranchType();
      createBranch(&defVal_);
    }
    virtual ~ExtendedBranch(){
      delete defVal;
    }

    // If theTree contains a branch br with name bname already, br->SetAddress(addr) needs to be called.
    // However, addr for a vector data has to be of type vector<varType>**, not vector<varType>*.
    // This necessitates making valueArray a pointer, which complicates the implementation considerably.
    // Therefore, copy and move constructors, or operators, are deleted.
    // This deletion also implies that one cannot call construct a vector<ExtendedBranch>.
    // A vector<ExtendedBranch*> is still permissible.
    // A default constructor does not make sense for this class either, so it is also deleted.
    ExtendedBranch() = delete;
    ExtendedBranch(ExtendedBranch const&) = delete;
    ExtendedBranch& operator=(ExtendedBranch const&) = delete;
    ExtendedBranch(ExtendedBranch&& other) = delete;
    ExtendedBranch& operator=(ExtendedBranch&& other) = delete;


  protected:

    void assignBranchType(){
      btype = hasType<varType>::type;
      if (btype == nBranchTypes){
        std::cerr << "ExtendedBranch::assignBranchType: Could not determine the branch type for branch " << bname << std::endl;
        assert(0);
      }
    }
    void createBranch(varType const* defVal_){
      if (!isVector && defVal_){
        delete defVal;
        defVal = new varType(*defVal_);
      }
      reset();

      TString leaftypename = "";
      if (!isVector){
        if (btype == BranchHelpers::bBool) leaftypename = "O";
        else if (btype == BranchHelpers::bUChar) leaftypename = "b";
        else if (btype == BranchHelpers::bChar) leaftypename = "B";
        else if (btype == BranchHelpers::bUShort) leaftypename = "s";
        else if (btype == BranchHelpers::bShort) leaftypename = "S";
        else if (btype == BranchHelpers::bUInt) leaftypename = "i";
        else if (btype == BranchHelpers::bInt) leaftypename = "I";
        else if (btype == BranchHelpers::bULongLong) leaftypename = "l";
        else if (btype == BranchHelpers::bLongLong) leaftypename = "L";
        else if (btype == BranchHelpers::bFloat) leaftypename = "F";
        else if (btype == BranchHelpers::bDouble) leaftypename = "D";
      }
      if (theTree && bname!=""){
        if (leaftypename==""){
          if (isVector) theTree->Branch(bname, &valueArray);
          else theTree->Branch(bname, &value);
        }
        else theTree->Branch(bname, &value, (bname + "/" + leaftypename));
      }
    }

  };

}

#endif
