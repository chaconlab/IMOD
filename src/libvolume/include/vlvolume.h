/***************************************************************************
                          vlvolume.h  -  description
                             -------------------
    begin                : Mon Apr 5 2004
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


#ifndef _vlVolume_h
#define _vlVolume_h


#include "vlenums.h"
#include "tripletypes.h"
//#include "vlvoliterbase.h"
#include "vlvoldatalayoutbase.h"


/**
 * THIS IS THE MAIN CLASS OF THE LIBRARY
 *
 * vlVolume class is the central class for a 3D volume. This class should be used
 * by users to work with 3D volume datasets. The class provides easy to use API
 * for reading/writing volume from/to disk, accessing volumetric data through the
 * use of iterators, accessing volume information etc.
 * This is the volume class of OpenVL visible to the users. This class brings together
 * all the different components : data layouts, iterators, file I/O filters etc.
 *
 * @author Sarang Lakare (lsarang@cs.sunysb.edu)
 */
class vlVolume {
// Member functions
public:

// all iterators should be able to access the internals of the volume class
  template <class T>
  friend class vlVolIterBase;

  template <class T>
  friend class vlVolIterBaseConst;

  /**
   *  default constructor
   *
   *  @param dim: Dimensions in voxels
   *  @param dataType: voxel type
   *  @param units: Voxel dimensions
   */
  vlVolume(const vlDim & dim = vlDim(1,1,1),
           const vlDataType dataType = UnsignedInt8,
           const vlUnit & units = vlUnit(1.0, 1.0, 1.0));

  /**
   * Constructor where the layout can be specified
  *
  *  @param dim: Dimensions in voxels
  *  @param layout: Layout type of the map
  *  @param dataType: voxel type
  *  @param units: Voxel dimensions
  */
  vlVolume(const vlDim & dim,
           const std::string & layout,
           const vlDataType dataType = UnsignedInt8,
            const vlUnit & units = vlUnit(1.0, 1.0, 1.0));

  /**
  * Copy Constructor
  * @param old: Reference Volume to be copied
  */
  vlVolume(vlVolume *old);


  /**
   *  Constructor to create a padded volume from a previous volume (padding is introduced in both sides of the volume)
   *  @param old: Reference Volume to be copied
   *  @param margin: Padded dimensions
   */
  vlVolume(vlVolume *old,vlDim margin);

   /**
    *  Constructor to create a incremented volume from a previous volume (padding is introduced in only on side of the volume)
    *  @param old: Reference Volume to be copied
    *  @param inc: Padded dimensions
    *  @param dummy: useless parameter
    */
  vlVolume(vlVolume *old,vlDim inc, int dummy);

  /// default destructor
  virtual ~vlVolume();

  ///Checks if this is a valid volume. Returns true of it is valid, else false.
  bool isValid() const;

  /// Clears the volume data with the given data
  bool clear(const uint8 data = 0x00);

  /// Returns the dimensions of the 3D data
  vlDim dim() const;

  /// Returns the stepping in voxel distance along x, y and z axis
  vlStep stepping() const;

  /// Returns the actual distance between the voxels along x, y and z
  vlUnit units() const;

  /// Returns the number of bits that form a voxel
  uint16 bitsPerVoxel() const;

  /// Returns the number of bytes that form a voxel
  uint16 bytesPerVoxel() const;

  /// Returns the total number of voxels in the volume = XDim x YDim x ZDim
  uint64 voxelCount() const;

  /// Returns the type of the data stored
  vlDataType dataType() const;

  /// Returns the type info of the stored data
  const char * typeInfo() const;

  /// Returns the layout in which data is stored
  vlLayoutType dataLayout() const;


  /// Returns the state of the dirty flag
  bool isDirty() const;

  /// Set/Reset the dirty flag
  void setDirty(bool dirty);

  /**
   * Returns the pointer to the voxel at the given position. This returns a void pointer,
   * so make sure you cast it to the correct type. To avoid type-conflicts, "use iterators".
   *
   * @param position: location of the voxel in 3D.
   * @return void * pointer to the memory location of the voxel.
   */
  void * getVoxelVoidPtr(const vlPoint3ui & position) const;

  /// Fast access to a voxel. This returns a void pointer,
  ///so make sure you cast it to the correct type.
  /// @param position:        location of the voxel in the list.
  float inline getVoxelSpecial(const int pos) const
  {
     return ((vlVolDataLayoutBase<float> *)m_pVolData)->getVoxelSpecial(pos);
  };

  ///Fast write of a voxel
  /// @param position:        location of the voxel in 3D.
  bool inline setVoxelSpecial(const vlPoint3ui & position, const float voxel)
  {
    ((vlVolDataLayoutBase<float> *)m_pVolData)->setVoxelSpecial(position, voxel);
    return (true);
  };

  ///Fast write of a voxel
  /// @param position:        location of the voxel in the list.
  /// @param voxel: 		  value to write
  template <class DataType>
  bool inline setVoxelSpecial(const int & position, const DataType voxel)
  {
    ((vlVolDataLayoutBase<DataType> *)m_pVolData)->setVoxelSpecial(position, voxel);
    return (true);
  };

  /// Get the voxel value at the position  of a voxel
  ///@param: Position of the voxel
  /// @param voxel: return the read value
  template <class DataType>
  bool getVoxel(const vlPoint3ui & position, DataType & voxel) const
  {
    if(m_pVolData) {
      if(typeid(DataType).name() == m_pVolData->typeInfo()) {
        voxel = ((vlVolDataLayoutBase<DataType> *)m_pVolData)->getVoxel(position);
        return (true);
      } else {
    	  fprintf(stderr,"ERROR : Invalid datatype in getVoxel(..) : Actual datatype = %s\n",m_pVolData->typeInfo());
        return (false);
      }
    } else
    return (false);
  };

  /// Get the voxel value at a real position. It makes an interpolation from the
  /// voxels around the position
  /// @param: Real position: Position of the voxel
  /// @param voxel: return the read value
  template <class DataType>
  bool getVoxel(const vlPoint3f & position, DataType & voxel) const
  {
    if(m_pVolData) {
      if(typeid(DataType).name() == m_pVolData->typeInfo()) {
        voxel = ((vlVolDataLayoutBase<DataType> *)m_pVolData)->getVoxel(position);
        return (true);
      } else {
    	  fprintf(stderr,"ERROR : Invalid datatype in getVoxel(..) : Actual datatype = %s\n",m_pVolData->typeInfo());
        return (false);
      }
    } else
      return (false);
  };

  /// Set the voxel value at the given position to 'voxel'
  template <class DataType>
  bool setVoxel(const vlPoint3ui & position, const DataType voxel)
	{
		if(m_pVolData) {
  		if(typeid(DataType).name() == m_pVolData->typeInfo()) {
	  	  /*voxel=*/((vlVolDataLayoutBase<DataType> *)m_pVolData)->setVoxel(position, voxel);
				return (true);
			} else {
				fprintf(stderr,"ERROR : Invalid datatype in getVoxel(..) : Actual datatype = %s\n",m_pVolData->typeInfo());
				return (false);
			}
		} else
			return (false);
	};

  /// Increase with an interpolation the voxels around the given position.
  /// ONLY IMPLEMENTED FOR FLOAT TYPE
  template <class DataType>
  bool interpolateVoxel(const vlPoint3f & position, const DataType voxel)
  {
                    return (true);
  };
  bool interpolateVoxel(const vlPoint3f & position, const float voxel)
  {
            ((vlVolDataLayoutBase<float> *)m_pVolData)->incVoxel(position, voxel);
            return (true);
  };



 /**
  * Returns the reference position of the volume in the space
 *
 *  @param: Tridimensional position returned
 */
  void getPosition(vlPoint3f *p);

  ///Set the reference position of the volumen in the space
  void setPosition(vlPoint3f p);

  /// Returns a const pointer to the storage of volume data.
  vlVolData * volData() const { return (m_pVolData); };

  ///Return the position in the space of the center of the volume
  vlPoint3f getCenter();

protected:

  /// Returns a RW pointer to the storage of volume data.
  vlVolData * volDataRW() { return (m_pVolData); };

private:

  /// Auxiliar Function to create the Volume data buffer
  template <class DataType>
  vlVolData * getLayoutT(DataType & dummy, const vlDim & dim, const vlUnit & units);

  /// stores the volume data
  vlVolData *m_pVolData;

  /// Volume position in the real space
  vlPoint3f position;
};

#endif // _vlVolume_h

