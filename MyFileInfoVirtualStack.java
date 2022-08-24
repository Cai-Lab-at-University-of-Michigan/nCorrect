//https://github.com/YuToyoshima/getSubImage/blob/master/MyFileInfoVirtualStack.java

import ij.*;
import ij.process.*;
import ij.io.*;
import java.io.*;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Properties;
import java.util.ArrayList;
import java.util.List;

/** This plugin opens a multi-page TIFF file as a virtual stack. It
 * implements the File/Import/TIFF Virtual Stack command. */
public class MyFileInfoVirtualStack extends VirtualStack {
    FileInfo[] info;
    int nImages;
    public RandomAccessFile raf;
    public List<ByteBuffer> mappings = new ArrayList<ByteBuffer>();
    public static final int MAPPING_SIZE_MAX = 1 << 30; // MappedByteBuffer cannot handle over 2GB
    public long mapping_size;
    public long sliceBytes;
    public long sliceBytesWithGap;
    ImagePlus myimp;
    
    
    /* Constructs a FileInfoVirtualStack from a FileInfo object. */
    public MyFileInfoVirtualStack(FileInfo fi) {
        info = new FileInfo[1];
        info[0] = fi;
        myOpen(true);
    }
    
    /* Constructs a FileInfoVirtualStack from a FileInfo
        object and displays it if 'show' is true. */
    public MyFileInfoVirtualStack(FileInfo fi, boolean show) {
        info = new FileInfo[1];
        info[0] = fi;
        myOpen(show);
    }
    
    public MyFileInfoVirtualStack(FileInfo[] info, boolean show) {
        this.info = info;
        myOpen(show);
    }
    
    public MyFileInfoVirtualStack(String path, boolean show) {
        File file = new File(path);
        try {
            raf = new RandomAccessFile(file,"r");
        } catch (FileNotFoundException e) {
            // raf.close();
            IJ.error("Virtual Stack","File Not Found");
        }
        
        TiffDecoder td = new TiffDecoder(file.getParent()+file.separator, file.getName());
        try {info = td.getTiffInfo();}
        catch (IOException e) {
            String msg = e.getMessage();
            if (msg==null||msg.equals("")) msg = ""+e;
            IJ.error("try-catch in TiffDecoder", msg);
            return;
        }
        if (info==null || info.length==0) {
            IJ.error("Virtual Stack", "This does not appear to be a TIFF stack");
            return;
        }
        if (IJ.debugMode)
            IJ.log(info[0].debugInfo);
        myOpen(show);
    }
    
    void myOpen(boolean show) {
        FileInfo fi = info[0];
        int n = fi.nImages;
        sliceBytes = fi.width*fi.height*fi.getBytesPerPixel();
        sliceBytesWithGap = sliceBytes + fi.gapBetweenImages;
        if (info.length==1 && n>1) {
            info = new FileInfo[n];
            
            for (int i=0; i<n; i++) {
                info[i] = (FileInfo)fi.clone();
                info[i].nImages = 1;
                info[i].longOffset = fi.getOffset() + i*sliceBytesWithGap;
            }
        }
        nImages = info.length;
        
        /* setup MappedByteBuffer */
        FileChannel fch = raf.getChannel();
        try {
            long fileSize = raf.length();
            mapping_size = (MAPPING_SIZE_MAX/sliceBytesWithGap)*sliceBytesWithGap;
            for (long offset=fi.getOffset(); offset<fileSize; offset+=mapping_size) {
                long size2 = Math.min(fileSize-offset, mapping_size);
                mappings.add(fch.map(FileChannel.MapMode.READ_ONLY,offset,size2));
            }} catch (IOException e) {
                IJ.error("Virtual Stack","Error occurred in setup MappedByteBuffer");
            }
        
        FileOpener fo = new FileOpener(info[0] );
        ImagePlus imp = fo.open(false);
        if (nImages==1 && fi.fileType==FileInfo.RGB48) {
            if (show) imp.show();
            myimp = imp;
            return;
        }
        Properties props = fo.decodeDescriptionString(fi);
        ImagePlus imp2 = new ImagePlus(fi.fileName, this);
        imp2.setFileInfo(fi);
        if (imp!=null && props!=null) {
            setBitDepth(imp.getBitDepth());
            imp2.setCalibration(imp.getCalibration());
            imp2.setOverlay(imp.getOverlay());
            if (fi.info!=null)
                imp2.setProperty("Info", fi.info);
            int channels = getInt(props,"channels");
            int slices = getInt(props,"slices");
            int frames = getInt(props,"frames");
            if (channels*slices*frames==nImages) {
                imp2.setDimensions(channels, slices, frames);
                if (getBoolean(props, "hyperstack"))
                    imp2.setOpenAsHyperStack(true);
            }
            if (channels>1 && fi.description!=null) {
                //int mode = IJ.COMPOSITE;
                int mode = CompositeImage.COMPOSITE;
                if (fi.description.indexOf("mode=color")!=-1)
                    //mode = IJ.COLOR;
                    mode = CompositeImage.COLOR;
                else if (fi.description.indexOf("mode=gray")!=-1)
                    //mode = IJ.GRAYSCALE;
                    mode = CompositeImage.GRAYSCALE;
                imp2 = new CompositeImage(imp2, mode);
            }
        }
        if (show) imp2.show();
        myimp = imp2;
    }
    
    
    public ImagePlus getImage() {
        return myimp;
    }
    
    
    int getInt(Properties props, String key) {
        Double n = getNumber(props, key);
        return n!=null?(int)n.doubleValue():1;
    }
    
    Double getNumber(Properties props, String key) {
        String s = props.getProperty(key);
        if (s!=null) {
            try {
                return Double.valueOf(s);
            } catch (NumberFormatException e) {}
        }
        return null;
    }
    
    boolean getBoolean(Properties props, String key) {
        String s = props.getProperty(key);
        return s!=null&&s.equals("true")?true:false;
    }
    
    /** Deletes the specified image, were 1<=n<=nImages. */
    public void deleteSlice(int n) {
        if (n<1 || n>nImages)
            throw new IllegalArgumentException("Argument out of range: "+n);
        if (nImages<1) return;
        for (int i=n; i<nImages; i++)
            info[i-1] = info[i];
        info[nImages-1] = null;
        nImages--;
    }
    
    /** Returns an ImageProcessor for the specified image,
     * were 1<=n<=nImages. Returns null if the stack is empty.
     */
    public ImageProcessor getProcessor(int n) {
        if (n<1 || n>nImages)
            throw new IllegalArgumentException("Argument out of range: "+n);
        if (IJ.debugMode) IJ.log("FileInfoVirtualStack: "+n+", "+info[n-1].getOffset());
        //if (n>1) IJ.log("  "+(info[n-1].getOffset()-info[n-2].getOffset()));
        info[n-1].nImages = 1; // why is this needed?
        FileOpener fo = new FileOpener(info[n-1]);
        ImagePlus imp = fo.open(false);
        if (imp!=null)
            return imp.getProcessor();
        else {
            int w=getWidth(), h=getHeight();
            IJ.log("Read error or file not found ("+n+"): "+info[n-1].directory+info[n-1].fileName);
            switch (getBitDepth()) {
                case 8: return new ByteProcessor(w, h);
                case 16: return new ShortProcessor(w, h);
                case 24: return new ColorProcessor(w, h);
                case 32: return new FloatProcessor(w, h);
                default: return null;
            }
        }
    }
    
    /** Returns the number of images in this stack. */
    public int getSize() {
        return nImages;
    }
    
    /** Returns the label of the Nth image. */
    public String getSliceLabel(int n) {
        if (n<1 || n>nImages)
            throw new IllegalArgumentException("Argument out of range: "+n);
        if (info[0].sliceLabels==null || info[0].sliceLabels.length!=nImages)
            return null;
        else
            return info[0].sliceLabels[n-1];
    }
    
    public int getWidth() {
        return info[0].width;
    }
    
    public int getHeight() {
        return info[0].height;
    }
       
    
    /* 
     * Return pixel values of subimage specified with slice and range.
     * This function works faster than getPixels by the following reasons.
     * 1. This function does not access unneccesary pixels.
     * 2. This function is specialized to get pixel values. 
     *    Unnecessary processes (i.e. LUT) are omitted.
     * 3. This function utilizes java.nio.MappedByteBuffer (MemoryMapped I/O).
     */
    public Object getSubImage(int slice, int[] range) {
        // public Object getPixelsPartial_mbb(int slice, int[] range) {
        int bytePerPixel = info[0].getBytesPerPixel();
        int width = info[0].width;
        int height = info[0].height;
        int offset_width = width*bytePerPixel;
        int x_start = range[0];
        int y_start = range[1];
        int x_end = range[2];
        int y_end = range[3];
        int width_subim  = x_end - x_start + 1;
        int height_subim = y_end - y_start + 1;
        int width_subim_byte = width_subim*bytePerPixel;
        
        long offset_to_image = info[slice-1].longOffset-info[0].getOffset();
        int idx_map       = (int) (offset_to_image / mapping_size);
        int offset_in_map = (int) (offset_to_image % mapping_size);
        MappedByteBuffer mbb = (MappedByteBuffer) mappings.get(idx_map);
        
        byte[] buf = new byte[width_subim_byte*height_subim];
        for (int cy=0; cy<=y_end-y_start; cy++) {
            mbb.position(offset_in_map + offset_width*(cy+y_start) + x_start*bytePerPixel);
            mbb.get(buf,width_subim_byte*cy,width_subim_byte);
        }
        
        switch (bytePerPixel) {
            case 1:
                return buf;
            case 2:
                return convertByte2Short(buf);
            case 4:
                return convertByte2Float(buf);
        }
        return null;
    }
    
    public Object getSubImage(int slice) {
        int width  = info[0].width;
        int height = info[0].height;
        int[] range = {0,0,width-1,height-1};
        return getSubImage(slice,range);
    }
    
/*    
    public Object getPixelsPartial_raf(int slice, int[] range) {
        int bytePerPixel = info[0].getBytesPerPixel();
        int width = info[0].width;
        int height = info[0].height;
        int offset_width = width*bytePerPixel;
        long offset_to_image = info[slice-1].longOffset;
        int x_start = range[0];
        int y_start = range[1];
        int x_end = range[2];
        int y_end = range[3];
        int width_subim  = x_end - x_start + 1;
        int height_subim = y_end - y_start + 1;
        int width_subim_byte = width_subim*bytePerPixel;
        
        byte[] buf = new byte[width_subim_byte*height_subim];
        try{
            for (int cy=0; cy<=y_end-y_start; cy++) {
                raf.seek(offset_to_image + offset_width*(cy+y_start) + x_start*bytePerPixel);
                raf.read(buf,width_subim_byte*cy,width_subim_byte);
            }
        } catch (IOException e) {
            IJ.error("Virtual Stack","Error occured in access RandomAccessFile");
        }
        
        switch (bytePerPixel) {
            case 1:
                return buf;
            case 2:
                return convertByte2Short(buf);
            case 4:
                return convertByte2Float(buf);
        }
        return null;
    }
    */
    
    /* Convert byte arrays to short arrays. Derived from ij.io.ImageReader.java */
    public Object convertByte2Short(byte[] buffer) {
        FileInfo fi = info[0];
        int base = 0;
        int bufferSize = buffer.length;
        int bytesPerPixel = fi.getBytesPerPixel();
        int pixelsRead = bufferSize/bytesPerPixel;
        short[] pixels = new short[pixelsRead];
        if (fi.intelByteOrder) {
            if (fi.fileType==FileInfo.GRAY16_SIGNED)
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)((((buffer[j+1]&0xff)<<8) | (buffer[j]&0xff))+32768);
            else
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)(((buffer[j+1]&0xff)<<8) | (buffer[j]&0xff));
        } else {
            if (fi.fileType==FileInfo.GRAY16_SIGNED)
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)((((buffer[j]&0xff)<<8) | (buffer[j+1]&0xff))+32768);
            else
                for (int i=base,j=0; i<(base+pixelsRead); i++,j+=2)
                    pixels[i] = (short)(((buffer[j]&0xff)<<8) | (buffer[j+1]&0xff));
        }
        return pixels;
    }
    
    /* Convert byte arrays to float arrays. Derived from ij.io.ImageReader.java */
    public Object convertByte2Float(byte[] buffer) {
        FileInfo fi = info[0];
        int base = 0;
        int bufferSize = buffer.length;
        int bytesPerPixel = fi.getBytesPerPixel();
        int pixelsRead = bufferSize/bytesPerPixel;
        float[] pixels = new float[pixelsRead];
        int pmax = base + pixelsRead;
        int tmp;
        int j = 0;
        if (fi.intelByteOrder)
            for (int i=base; i<pmax; i++) {
            tmp = (((buffer[j+3]&0xff)<<24) | ((buffer[j+2]&0xff)<<16) | ((buffer[j+1]&0xff)<<8) | (buffer[j]&0xff));
            if (fi.fileType==FileInfo.GRAY32_FLOAT)
                pixels[i] = Float.intBitsToFloat(tmp);
            else if (fi.fileType==FileInfo.GRAY32_UNSIGNED)
                pixels[i] = (float)(tmp&0xffffffffL);
            else
                pixels[i] = tmp;
            j += 4;
            }
        else
            for (int i=base; i<pmax; i++) {
            tmp = (((buffer[j]&0xff)<<24) | ((buffer[j+1]&0xff)<<16) | ((buffer[j+2]&0xff)<<8) | (buffer[j+3]&0xff));
            if (fi.fileType==FileInfo.GRAY32_FLOAT)
                pixels[i] = Float.intBitsToFloat(tmp);
            else if (fi.fileType==FileInfo.GRAY32_UNSIGNED)
                pixels[i] = (float)(tmp&0xffffffffL);
            else
                pixels[i] = tmp;
            j += 4;
            }
        return pixels;
    }
    
}