/***************************************************************************
                          vlvoliterbase.h  -  description
                             -------------------
    begin                : Mon Apr 12 2004
    copyright            : (C) 2004 by Jose Ignacio Garzon
    email                :
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/


 #ifndef _vlVolIterBase_h
#define _vlVolIterBase_h

#include <typeinfo>

#include "vlenums.h"
#include "tripletypes.h"
#include "vlvoldatalayoutbase.h"

template <typename DataType>
class vlInterpolatorBase;

template <typename DataType>
class vlVoxelOpBase;

class vlVolume;
class vlVoxelOpValue;


class vlOffset;

class vlNeighborhood;


/**
 * This is base class of all iterators - hence the name SuperBase. This class
 * is pretty much insignificant from a user's point of view. It is needed
 * to pass base pointer to iterator objects. For e.g., this is used in data layout
 * plugins to get pointer to the iterator object.
 */
class vlVolIterSuperBase {
public:
	/// Returns true if the iterator was properly initialized
	virtual bool isValid() = 0;
};


/**
 * This is the base class for volume iterators. Each layout type should
 * re-implement all the API functions in the most optimal way.
 * Iterator class for each data layout should have the name vlVolIter<layout_name>.
 * For e.g., for Linear layout, the iterators are vlVolIterLinear and vlVolIterLinearConst.
 * This base class if for const iterators. That is, iterators which will not modify
 * any data. You will also need to implement iterator derived from the non-const version
 * which is vlVolIterBse.
 *
 * @see vlVolIterBase
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 */
template <typename DataType>
class vlVolIterBaseConst : public vlVolIterSuperBase {

public:
	/// Default constructor
	vlVolIterBaseConst() { m_volData = 0L; };

	/// Constructor when a const volume pointer is given
	vlVolIterBaseConst(const vlVolume * vol)
	{
		if(vol->typeInfo() == typeid(DataType).name()) {
			m_volData = dynamic_cast< vlVolDataLayoutBase<DataType>* >(vol->volData());
		} else {
			fprintf(stderr,"ERROR :: vlVolIterBaseConst<DataType> : Trying to call iterator with wrong datatype.\n");
			m_volData = 0L;
		}
	};

	/// Constructor when a const volume data pointer is given
	vlVolIterBaseConst(const vlVolData * data)
	{
		if(data->typeInfo() == typeid(DataType).name()) {
			const vlVolDataLayoutBase<DataType>* volData = dynamic_cast< const vlVolDataLayoutBase<DataType>* >(data);
			m_volData = const_cast< vlVolDataLayoutBase<DataType>* > (volData);
		} else {
			fprintf(stderr,"ERROR :: vlVolIterBaseConst<DataType> : Trying to call iterator with wrong datatype.\n");
			m_volData = 0L;
		}
	};

	/// Constructor when a const data buffer pointer is given
	vlVolIterBaseConst(const vlVolDataLayoutBase<DataType> * data)
	{
		m_volData = const_cast< vlVolDataLayoutBase<DataType>* > (data);
	};

	/// Default destructor
	virtual ~vlVolIterBaseConst() { };

  ///Return the volume object
  virtual vlVolDataLayoutBase<DataType> *getVolume() const=0;

  /// move the iterator to the first voxel
  virtual bool begin() = 0;

  /// Return true if at the end of the volume
  virtual bool end() = 0;

  /// Return the position of the iterator
  virtual vlPoint3ui pos() = 0;

  /// move the iterator to voxel at newPosition - for random access
  virtual bool moveTo(const vlPoint3ui & newPosition) = 0;

  /// move the iterator by offset
  virtual bool moveRelative(const vlOffset & offset) = 0;

  /// Get the value of voxel at current position
  virtual DataType get() = 0;

  /// Get the value of voxel at the given position
  virtual DataType getValueAt(const vlPoint3ui & pos) = 0;

  /**
   * Get the value of voxel at the given position - interpolation will be used.
   * Set the type of interpolation to use using setInterpolation().
   */
  virtual DataType getValueAt(const vlPoint3f & pos) = 0;

  /// move the iterator to the next position - same as next()
  virtual bool operator++() = 0;

  /// move the iterator to the previous position - same as prev()
  virtual bool operator--() = 0;

  /// move the iterator to the next voxel - NOTE : No order assumed - fastest traversal
  virtual bool next() = 0;

  /**
   * Same as next() but w/o checking for any bounds.
   * WARNING : Using this w/o proper thought is dangerous. You have been warned!
   * This is simply the fastest way to move to the next voxel. You need to do
   * bounds-checking yourself, else the pointer can easily go outside the data range.
   *
   * @see next()
   */
  virtual void nextNBC() = 0;

  /// move the iterator to the next voxel - move first along X, then Y and then Z
  virtual bool nextXYZ() = 0;

  /// move the iterator to the next voxel - move first along Y, then Z and then X
  virtual bool nextYZX() = 0;

  /// move the iterator to the next voxel - move first along Z, then X and then Y
  virtual bool nextZXY() = 0;

  /// move the iterator to the previous voxel - NOTE : No order assumed - fastest traversal
  virtual bool prev() = 0;

  /**
   * Same as prev() but w/o checking for any bounds.
   * WARNING : Using this w/o proper thought is dangerous. You have been warned!
   * This is simply the fastest way to move to the previous voxel. You need to do
   * bounds-checking yourself, else the pointer can easily go outside the data range.
   *
   * @see prev()
   */
  virtual void prevNBC() = 0;

  /// move the iterator to the previous voxel - move first along X, then Y and then Z
  virtual bool prevXYZ() = 0;

  /// move the iterator to the previous voxel - move first along Y, then Z and then X
  virtual bool prevYZX() = 0;

  /// move the iterator to the previous voxel - move first along Z, then X and then Y
  virtual bool prevZXY() = 0;

  /// get value at a voxel with position delta from the current one
  virtual DataType getRelative(const vlOffset & offset) = 0;

  /// Same as getRelative but with NBC = No Boundary Check, thus faster.
  virtual DataType getRelativeNBC(const vlOffset & offset) = 0;

  /**
  * Get value at a point delta from the curret position of the iterator.
  * This will be an iterpolated value. Set the type of interpolation
  * using setIterpolation()
  */
  virtual DataType getRelative(const vlPoint3f & delta) = 0;

  /// Get the voxel along X axis
  virtual DataType getRelativeX(const int32 offset) = 0;

  /// Get the voxel along Y axis
  virtual DataType getRelativeY(const int32 offset) = 0;

  /// Get the voxel along Z axis
  virtual DataType getRelativeZ(const int32 offset) = 0;

  /**
  * Returns true if the outside flag was set. The outside flag is set whenever
  * a get*() function is called (except the NBC versions). It is set to true if
  * the get*() tried to access outside the volume (in which case it returned zero).
  * It is set to false, if the get call was successful. Thus, by querying this flag,
  * you can find out if get*() was inside volume. The reason why outside flag works
  * only for get*() functions is that these functions have no way to return if
  * they were successful or not (since they return voxel value). All other functions
  * return bool indicating if the function was successful.
  */
  virtual bool outside() = 0;

  /// Explicitly reset the outside flag
  virtual void resetOutside() = 0;

  /**
   * Set the neighborhood of the current voxel to access. A neighborhood
   * is defined by a list of offsets stored in the vlNeighborhood object.
   *
   * @see vlNeighborhood
   */
  virtual bool setNeighborhood(const vlNeighborhood & neighborhood) = 0;

  /**
   * Get the value stored at currently selected neighboring voxel. If a
   * neighborhood is setup using a call to setNeighborhood(..), then
   * there will always be a neighboring voxel that is selected at any given
   * time. This function call will return the value at that voxel.
   */
  virtual DataType getNeighbor() = 0;

  /**
   * Get the offset of the current neighbor in the neighborhood. The offset
   * is the location of the neighbor w.r.t. the current location.
   */
  virtual vlOffset getNeighborOffset() = 0;

  /**
   * Get the id of the current neighbor in the neighborhood. This id is basically
   * the position of the neighbor in the neighborhood list. Thus, this id can be
   * used to access the neighbor's properties stored in vlNeighborhood or the classes
   * that derive from it.
   */
  virtual uint16 getNeighborId() = 0;

  /**
   * Move the neighbor pointer to the next neighbor in the neighborhood list. If
   * there are no more neighbors left, then the function will return false.
   */
  virtual bool nextNeighbor() = 0;

  /**
   * Move the neighbor pointer to the first neighbor in the neighborhood list. Call
   * this function to restart looping over the neighboring voxels.
   */
  virtual bool firstNeighbor() = 0;

  /**
   * This function checks if end of neighborhood is reached. If it has reached, it
   * returns true, else it returns false. End of neighborhood is defined as a call
   * to nextNeighbor() when there are no more neighbors. That is, when nextNeighbor()
   * is called while at the last neighbor. This function is useful mainly to write
   * neat and easy-to-understand code like this:
   * @code
   * while(!iter.endOfNeighborhood()) {
   *  std::cout << iter.getNeighbor() << std::flush;
   *  iter.nextNeighbor();
   * }
   * @endcode
   */
  virtual bool endOfNeighborhood() = 0;

  /***
  * Set the Interpolator that will be used to
  * compute values in real positions
  */
  virtual bool setInterpolation(vlInterpolatorBase<DataType> *interP) = 0;

  /**
   * Returns pointer to the interpolator being used. Returns 0L if
   * no interpolator is set. This pointer can be used to get the
   * type and name of the interpolator.
   */
   virtual vlInterpolatorBase<DataType> * getInterpolator() = 0;

  /**
  * Sets the Operator of the Iterator
  */
  virtual bool setVoxelOperation(vlVoxelOpBase<DataType> *op) = 0;

  /**
  * Returns a pointer to the operator of the Iterator
  */
  virtual vlVoxelOpBase<DataType> const * getVoxelOperator() = 0;

  /**
  * Returns the value given by the operator
  * */
  virtual bool getOp(vlVoxelOpValue & value) = 0;

  /**
  * Writes the Volume in a File
  * */
  // virtual bool write(char *fileName,vlPoint3f *position)= 0 ;


protected:
	/// Pointer to the volumetric data of the volume
	vlVolDataLayoutBase<DataType> *m_volData;

};











/**
 * This is the base class for non-const volume iterators. Each layout type should
 * re-implement all the API functions in the most optimal way. Derive non-const
 * iterator for your data layout from this class. You will also need to implement
 * const iterators. Those will be derived from vlVolIterBaseConst.
 * Establece las declaraciones virtuales para un iterador de los metodos
 * QUE MODIFICAN los datos de un volumen
 * @see vlVolIterBaseConst
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 */
template <class DataType>
class vlVolIterBase { //: public vlVolIterBaseConst<DataType> {
public:
  /// Set the value of voxel at current position
  virtual bool set(const DataType & value) = 0;

  /// set value at a voxel offset from the current one
  virtual bool setRelative(const vlOffset & offset, const DataType & value) = 0;

	/// Same as getRelative but with NBC = No Boundary Check, thus faster.
	virtual bool setRelativeNBC(const vlOffset & offset, const DataType & value) = 0;

	/// Set the voxel along X axis
	virtual bool setRelativeX(const int32 offset, const DataType & value) = 0;

	/// Set the voxel along Y axis
	virtual bool setRelativeY(const int32 offset, const DataType & value) = 0;

	/// Set the voxel along Z axis
	virtual bool setRelativeZ(const int32 offset, const DataType & value) = 0;

  /// Set the current neighbor of the current voxel.
  virtual bool setNeighbor(const DataType & value) = 0;

}; // END - class VolIterBase









/**
 * This is the base class for non-const volume iterators. Each layout type should
 * re-implement all the API functions in the most optimal way. Derive non-const
 * iterator for your data layout from this class. You will also need to implement
 * const iterators. Those will be derived from vlVolIterBaseConst.
 * Establece las declaraciones virtuales para un iterador de los metodos
 * NO QUE MODIFICAN los datos de un volumen. Solo establece los metodos que no modifican
 * el volumen
 * Es un template para las clases iteradores donde se debera de establecer
 * el tipo de dato del volumen y el layout de este
 * @see vlVolIterBaseConst
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 */

template <typename DataType, vlLayoutType Layout = vlLayout::DefaultLayout>
class vlVolIterConst : public vlVolIterBaseConst<DataType> {

public:
	/// Default constructor
	vlVolIterConst(const vlVolume * vol)
		: vlVolIterBaseConst<DataType>(vol)
	{
	};


	vlVolIterConst(const vlVolData * data)
		: vlVolIterBaseConst<DataType>(data)
	{
	};


	/// Default destructor
	virtual ~vlVolIterConst() { };

  vlVolDataLayoutBase<DataType> *getVolume() const;

 /// move the iterator to the first voxel
  bool begin();

  /// Return true if at the end of the volume
  bool end() ;

  /// Return the position of the iterator
  vlPoint3ui pos();

  /// move the iterator to voxel at newPosition - for random access
  bool moveTo(const vlPoint3ui & newPosition);

  /// move the iterator by offset
  bool moveRelative(const vlOffset & offset);

  /// Get the value of voxel at current position
  DataType get();

  /// Get the value of voxel at the given position
  DataType getValueAt(const vlPoint3ui & pos);

  /**
   * Get the value of voxel at the given position - interpolation will be used.
   * Set the type of interpolation to use using setInterpolation().
   */
  DataType getValueAt(const vlPoint3f & pos);

  /// move the iterator to the next position - same as next()
  bool operator++();

  /// move the iterator to the previous position - same as prev()
  bool operator--();

  /// move the iterator to the next voxel - NOTE : No order assumed - fastest traversal
  bool next();

  /**
   * Same as next() but w/o checking for any bounds.
   * WARNING : Using this w/o proper thought is dangerous. You have been warned!
   * This is simply the fastest way to move to the next voxel. You need to do
   * bounds-checking yourself, else the pointer can easily go outside the data range.
   *
   * @see next()
   */
  void nextNBC();

  /// move the iterator to the next voxel - move first along X, then Y and then Z
  bool nextXYZ();

  /// move the iterator to the next voxel - move first along Y, then Z and then X
  bool nextYZX();

  /// move the iterator to the next voxel - move first along Z, then X and then Y
  bool nextZXY();

  /// move the iterator to the previous voxel - NOTE : No order assumed - fastest traversal
  bool prev();

  /**
   * Same as prev() but w/o checking for any bounds.
   * WARNING : Using this w/o proper thought is dangerous. You have been warned!
   * This is simply the fastest way to move to the previous voxel. You need to do
   * bounds-checking yourself, else the pointer can easily go outside the data range.
   *
   * @see prev()
   */
  void prevNBC();

  /// move the iterator to the previous voxel - move first along X, then Y and then Z
  bool prevXYZ();

  /// move the iterator to the previous voxel - move first along Y, then Z and then X
  bool prevYZX();

  /// move the iterator to the previous voxel - move first along Z, then X and then Y
  bool prevZXY();

  /// get value at a voxel with position delta from the current one
  DataType getRelative(const vlOffset & offset);

	/// Same as getRelative but with NBC = No Boundary Check, thus faster.
	DataType getRelativeNBC(const vlTriple<int32> & delta);

	/**
	 * Get value at a point delta from the curret position of the iterator.
	 * This will be an iterpolated value. Set the type of interpolation
   * using setIterpolation()
	 */
	DataType getRelative(const vlPoint3f & delta);

	/// Get the voxel along X axis
	DataType getRelativeX(const int32 offset);

	/// Get the voxel along Y axis
	DataType getRelativeY(const int32 offset);

	/// Get the voxel along Z axis
	DataType getRelativeZ(const int32 offset);

  bool setNeighborhood(const vlNeighborhood & neighborhood);

  DataType getNeighbor();

  vlOffset getNeighborOffset();

  uint16 getNeighborId();

  bool nextNeighbor();

  bool firstNeighbor();

  bool endOfNeighborhood();

	bool outside();

	void resetOutside();

 //Escribe volumen en fichero
  bool write(char *fileName,vlPoint3f *position) ;

  bool getOp(vlVoxelOpValue & value);
};







/**
 * This is the base class for non-const volume iterators. Each layout type should
 * re-implement all the API functions in the most optimal way. Derive non-const
 * iterator for your data layout from this class. You will also need to implement
 * const iterators. Those will be derived from vlVolIterBaseConst.
 * Esta clase combina las clases iteradoras de modificacion y no modificacion de
 * los datos de volumen
 * Es un template para las clases iteradores donde se debera de establecer
 * el tipo de dato del volumen y el layout de este
 *
 * @see vlVolIterBaseConst
 *
 * @author Sarang Lakare <sarang@users.sf.net>
 */
template <typename DataType, vlLayoutType Layout = vlLayout::DefaultLayout>
class vlVolIter : public vlVolIterConst<DataType, Layout>,
									 public vlVolIterBase<DataType> {

public:
	/// Default constructor
	vlVolIter(vlVolume * vol)
		: vlVolIterConst<DataType, Layout>(vol)
	{
	};

	vlVolIter(vlVolData * data)
		: vlVolIterConst<DataType, Layout>(data)
	{
	};

	/// Default destructor
	virtual ~vlVolIter() { };



  /// Set the value of voxel at current position
  bool set(const DataType & value);

  /// set value at a voxel offset from the current one
  bool setRelative(const vlOffset & offset, const DataType & value);

	/// Same as getRelative but with NBC = No Boundary Check, thus faster.
	bool setRelativeNBC(const vlOffset & offset, const DataType & value);

	/// Get the voxel along X axis
	bool setRelativeX(const int32 offset, const DataType & value);

	/// Get the voxel along Y axis
	bool setRelativeY(const int32 offset, const DataType & value);

	/// Get the voxel along Z axis
	bool setRelativeZ(const int32 offset, const DataType & value);

  bool setNeighbor(const DataType & value);

}; // END - class VolIterBase


#endif // _vlVolIterBase_h
