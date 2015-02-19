"""
 Script para interpolar los datos de mercator a la malla de Nemo
 usando la herramienta sosie.x para la interpolacion 3D.

 By Favio Medrano 
  2014-08-18 

"""

import os 
import sys 
import logging as log 
import datetime as dt 
import netCDF4 as nc 
import numpy as np
from os.path import join
# My imports
from netcdfFile import netcdfFile

SOSIE_EXE = 'sosie.x'

# Clase util

class fileEdit:
    """
    """
    handler = None
    content = None
    filename = None
    def __init__(self,fname):
        """
         Constructor de la clase, que lee el archivo de texto "fname"
         a modificar.
        """
        try:
            self.handler = open(fname,'r')
            self.filename = fname
        except Exception:
            log.exception('Fallo al abrir ' + fname)
            return None
        self.content = self.handler.readlines()
        self.handler.close()

    def ins(self,tag,nstr,jumpline=True):
        """
         Inserta una cadena delante en donde encuentre la cadena "tag"
        """
        tmp = []
        strj = '\n' if jumpline else ''
        for line in self.content:
            tmp.append(line)
            if tag in line:
                tmp.append(nstr + strj)
        self.content = tmp

    def rep(self,tag,nstr):
        """
         Reemplaza la(s) ocurrencia(s) de tag en el archivo por nstr
        """
        tmp = []
        for line in self.content:            
            if tag in line:
                tmp.append(line.replace(tag,nstr))
            else:
                tmp.append(line)
        self.content = tmp

    def saveFile(self,newfile=None):
        """
         Salva los cambios al archivo especificado en newfile, o hace una copia
         del archivo original (filename+'~') y salva el contenido en "filename"
        """
        if newfile == None:
            shutil.move(self.filename,self.filename+'~')
            self.handler = open(self.filename,'w')
        else:
            self.handler = open(newfile,'w')
        self.handler.writelines(self.content)
        self.handler.close()

#
# fileEdit class EOF
# --------------------------------------------------------

def set_unlimited(sRutaMercator):
    if not os.path.exists(sRutaMercator):
        log.error('No existe el archivo ' + sRutaMercator)
        return None

    dst = nc.Dataset(sRutaMercator,'r')
    sRutaMercatorTmp = sRutaMercator + 'tmp.nc'
    if not dst.dimensions['time_counter'].isunlimited():
        # Convertir la dimension time_counter sRutaMercator a unlimited
        # Usando ncks de nco 
        sCmdSetU = 'ncks -h --mk_rec_dmn time_counter ' + sRutaMercator + ' -o ' + sRutaMercatorTmp
        log.info('Comando: ' + sCmdSetU)
        os.system(sCmdSetU)

        if os.path.exists(sRutaMercatorTmp):
            log.info('Se genero el archivo ' + sRutaMercatorTmp)
            os.system('mv ' + sRutaMercatorTmp + ' ' + sRutaMercator)            
        else:
            log.error('Fallo el comando: ' + sCmdSetU)
            return None
    else:
        log.info('La dimension time_counter ya es unlimited')

    dst.close()
    return sRutaMercator


def total_seconds(tDelta):
    return abs(tDelta.seconds) + abs(tDelta.days) * 24 * 3600 


def extract_record(sRutaMercator,tFrame):
    if not os.path.exists(sRutaMercator):
        log.error('No existe el archivo ' + sRutaMercator)
        return None

    dst = nc.Dataset(sRutaMercator,'r')
    dates = nc.num2date(dst.variables['time_counter'][:], units = dst.variables['time_counter'].units, calendar = dst.variables['time_counter'].calendar) 
    ts = [ total_seconds(dates[i] - tFrame) for i in range(len(dates)) ]
    idx = np.argmin(ts)
    log.info('Buscando fecha : ' + str(tFrame)) 
    log.info('Indice a extraer ' + str(idx) + ' con fecha : ' + str(dates[idx])) 
    dst.close()
    sOutNc = 'mercator_data_' + dates[idx].strftime('%Y%m%d') + '.nc'
    sCmd = 'ncks -O -d time_counter,' + str(idx) + ',' + str(idx) + ' ' + sRutaMercator + ' -o ' + sOutNc 
    log.info('Comando para extraccion : ' + sCmd) 
    os.system(sCmd) 

    if os.path.exists(sOutNc):
        log.info('Se genero el archivo : ' + sOutNc)
        return sOutNc
    else:
        log.error('Fallo el comando ' + sCmd )
        return None 


def do_sosie_interp(sRutaMercator, sRutaMascara, sDate):

    # Verficar que existen los archivos
    if not (os.path.exists(sRutaMercator) and os.path.exists(sRutaMascara)):
        log.error('Alguno de los archivos Mercator o Mascara no existen, verificar.')
        return -1

    log.info('Datos mercator: ' + sRutaMercator)
    log.info('Archivo mascara: ' + sRutaMascara)    

    set_unlimited(sRutaMercator)
    # Extraer solo el registro que interesa
    outNc = extract_record(sRutaMercator,sDate)


    # Modificar los archivos namelist para usar sRutaMercator y sRutaMascara 
    namelistFiles = [ 'namelist.myocean.mercator.temp' , 'namelist.myocean.mercator.salinity' , \
                      'namelist.myocean.mercator.ssh' , 'namelist.myocean.mercator.u' , \
                      'namelist.myocean.mercator.v' ]

    outfiles = [ 'votemper_mercator-nemomesh_interpol.nc' , 'vosaline_mercator-nemomesh_interpol.nc' , \
                 'sossheig_mercator-nemomesh_interpol.nc' , 'vozocrtx_mercator-nemomesh_interpol.nc' , \
                 'vomecrty_mercator-nemomesh_interpol.nc' ]

    for idx,fNamelist in enumerate(namelistFiles):
        if os.path.exists(outfiles[idx]):
            log.info('El archivo ' + outfiles[idx] + ' ya existe, no se generara de nuevo. Skipping.')
            continue
        fModName = fileEdit(fNamelist)
        fModName.rep('INPUTFILE',outNc) 
        fModName.rep('MASKFILE',sRutaMascara) 
        fModName.rep('IDX','0')
        fModName.saveFile('namelist.ready')

        # Do sosie.x
        sCmd = SOSIE_EXE + ' -f namelist.ready' 
        log.info ('SOSIE CMD: ' + sCmd)
        os.system(sCmd)
        os.remove('namelist.ready') 
    #
    # Remove the extracted time netcdf
    os.remove(outNc)

    for f in outfiles:
        if not os.path.exists(f):
            log.error('Fallo la creacion del archivo: ' + str(f))
            return None 

    return outfiles


def make_bckIncFile(InData):

    t,z,y,x = InData['bckint'].shape
    # 
    # Definition of bck_assimilation_increments file 
    bckFileName = 'assim_background_increments.nc' 
    defDims = {'x': x, 'y': y , 't' : None , 'z' : z} 
    defVars = {'time_counter' : {'dimensions': ['t'],'dataType':'f8'} , 'time' : {'dimensions': [],'dataType':'f8'} , \
               'z_inc_dateb' : {'dimensions': [],'dataType':'f8'} , 'z_inc_datef' : {'dimensions': [],'dataType':'f8' } , \
               'bckint' : {'dimensions':['t','z','y','x'],'dataType':'f4'} , 'bckins' : {'dimensions':['t','z','y','x'],'dataType':'f4'} , \
               'bckineta' : {'dimensions':['t','y','x'],'dataType':'f4'} , 'bckinv' : {'dimensions':['t','z','y','x'],'dataType':'f4'} , \
               'bckinu' : {'dimensions':['t','z','y','x'],'dataType':'f4'}}

    bckFile = netcdfFile() 
    bckFile.createFile(bckFileName,'./','NETCDF3_CLASSIC')
    bckFile.createDims(defDims)
    bckFile.createVars(defVars) 
    bckFile.saveData(InData) 
    # TODO : Que va en las siguientes variables ??
    bckFile.fileHandler.variables['time'][:] = np.array([0])
    bckFile.fileHandler.variables['time_counter'][:] = np.array([0])
    bckFile.fileHandler.variables['z_inc_dateb'][:] = np.array([1])
    bckFile.fileHandler.variables['z_inc_datef'][:] = np.array([6])
    #
    bckFile.closeFile()
    return bckFileName



def do_bckIncData(filesDataMer, fileT, fileU, fileV, fDate = None):

    # Encontrar el indice mas cercano a fDate.
    if (fDate == None):
         timeF = 0
    else:
        dst = nc.Dataset(fileT,'r')
        log.info('Buscando indice cercano a : ' +str(fDate))
        timeF = nc.date2index(fDate,dst.variables['time_average_1d'], calendar=dst.variables['time_average_1d'].calendar, select='after')
        dst.close()
    log.info('fileT, fileU, fileV se lee el indice : ' + str(timeF))
    # Hacer la diferencia entre filesDataMer y filesModel para generar
    # el archivo bck_assimilation_increments.nc 

    varsInfo = ['votemper','vosaline','sossheig','vozocrtx','vomecrty'] 
    XAdataDic = {}
    XFdataDic = {}
    bckInc = {}
    # Load variables from XA and XF data
    for vName in varsInfo:
        #
        # XA Data
        for fileName in filesDataMer:
            if vName in fileName:
                dst = nc.Dataset(fileName,'r')
                XAdataDic[vName] = dst.variables[vName][:]
                dst.close()
        # 
        # XF Data
        if vName == 'votemper' or vName == 'vosaline' or vName == 'sossheig':
            fileToOpen = fileT 
        if vName == 'vozocrtx':
            fileToOpen = fileU
        if vName == 'vomecrty':
            fileToOpen = fileV
        dst = nc.Dataset(fileToOpen,'r')
        XFdataDic[vName] = dst.variables[vName][timeF]
        dst.close()
    # 
    # Convertir XA Data votemper a Celsius  
    XAdataDic['votemper'][ XAdataDic['votemper'][:]!=0 ] = XAdataDic['votemper'][ XAdataDic['votemper'][:]!=0 ] - 272.15 

    #
    # Hacer la diferencia. bckInc = XA - XF
    # 
    bckInc['bckint'] = XAdataDic['votemper'] - XFdataDic['votemper']
    bckInc['bckins'] = XAdataDic['vosaline'] - XFdataDic['vosaline']
    bckInc['bckint'][0,70:,:,:] = 0 
    bckInc['bckins'][0,70:,:,:] = 0
    bckInc['bckinu'] = (XAdataDic['vozocrtx'] - XFdataDic['vozocrtx']) * 0.8
    bckInc['bckinv'] = (XAdataDic['vomecrty'] - XFdataDic['vomecrty']) * 0.8
    bckInc['bckinu'][0,70:,:,:] = 0 
    bckInc['bckinv'][0,70:,:,:] = 0   
    bckInc['bckineta'] = XAdataDic['sossheig'] - XFdataDic['sossheig']
    #
    # Generar el archivo bck_assimilation_increments.nc 
    outF = make_bckIncFile(bckInc) 
    #
    if os.path.exists(outF):
        log.info('Se genero el archivo : ' + outF) 
        # 
        # Borrar archivos filesDataMer, pues son creados por sosie 
        # exclusivamente para esta operacion.
        for f in filesDataMer:
            os.remove(f) 
    else:
        log.error('Fallo al generarse el archivo : ' + outF) 

    return 1



def main():
    log.getLogger().setLevel(10)
    # Arguments en linea de comandos: 
    #  ruta a datos mercator, ruta archivo mascara

    if len(sys.argv) < 3:
        log.error('Missing Arguments') 
        log.error('\nUsage:\n\n > python make_iau_file.py $InputData $MaskFile\n')
        log.error('Loading defaults..')
        #return -1
        sRutaMercator = '/LUSTRE/hmedrano/STOCK/FORCING-RAW/MERCATOR_RAW/20140907/downloaded_mercator_query_20140826-20140914.nc'
        sRutaMascara = '/LUSTRE/hmedrano/STOCK/NEMO3.6/GOLFO12/GOLFO12-I/mesh_mask.nc'
    else:    
        sRutaMercator = sys.argv[1]
        sRutaMascara = sys.argv[2] 

    # Interpolar datos de entrada a la malla definida en sRutaMascara
    sosieFiles = do_sosie_interp(sRutaMercator,sRutaMascara,dt.datetime(2014,9,6))
    # sosieFiles contiene las variables interpoladas, votemper, vosaline, sossheig, vozocrtx, vomecrty
    #

    # Creamos el archivo assim_background_increments.nc , haciendo la diferencia entre los datos
    # de sosieFiles menos los del modelo.
    pathmodeloutputs = '/LUSTRE/hmedrano/STOCK/NEMO3.6/GOLFO12/GOLFO12-EXP001-S'
    tfile = 'GOLFO12-EXP001_1d_20140901_20140906_grid_T.nc'
    ufile = 'GOLFO12-EXP001_1d_20140901_20140906_grid_U.nc'
    vfile = 'GOLFO12-EXP001_1d_20140901_20140906_grid_V.nc'
    do_bckIncData(sosieFiles, join(pathmodeloutputs,tfile) , \
                              join(pathmodeloutputs,ufile) , \
                              join(pathmodeloutputs,vfile), dt.datetime(2014,9,6))
    # De salida tenemos el archivo: assim_background_increments.nc







if __name__ == "__main__":
    main()
