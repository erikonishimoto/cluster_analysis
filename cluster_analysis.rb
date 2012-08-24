require 'narray'

#-- calculating the value of the sum-of-square error
#--   for every possible merger of two clusters
def euclidean_distance( data, dim=0 )

  if dim == 1 or dim == -1
    data = data.transpose( dim, 0 )
  end

  nx,nt = data.shape

  imax = (nx-1)*nx/2  # 組み合せの数 nP2
  deu = NArray.float(imax)
  class_comb = Hash::new
  comb_class = Hash::new
  (0..nx-2).each{|i|
    (i+1..nx-1).each{|j|
      l = nx*i - (i*(i+1))/2 + j - i - 1
      d = data[i,true] - data[j,true]
      deu[l] = d.mul_add(d,0)

      #p [i,j,l]
      class_comb[l] = [i,j]
      comb_class[[i,j]] = l
    }
  }

  return deu,class_comb,comb_class
end


#-- Cluster Analysis with Ward's method
def cluster_analysis( data, dim, nclass=10, renumber=false )

  if dim == 1 or dim == -1
    data = data.transpose( dim, 0 )
  end
  nx,nt = data.shape

  #-- 各格子点の個体情報の初期化。はじめは全て違うクラス
  element_class = NArray.int(nx).indgen

  element_class_ary = NArray.int(nx,nx-nclass+1)
  element_class_ary[true,0] = element_class

  #-- 各格子間の平方ユークリッド距離の計算
  deu, class_comb, comb_class = euclidean_distance( data, 0 )

  #-- クラスタリング: 最終的にnclass個まで融合
  dmin_list = NArray.float(nx-nclass)
  nstep=0
  (0..nx-nclass-1).each{|istep|

    #-- 距離の最小値の算定
    dmin = deu.min
    dmin_list[nstep] = dmin
    #-- どのクラスの組み合わせか
    l = deu.eq(dmin).where[0]
    h,k = class_comb[l] #--> これはクラス番号

    #-- 新たに融合されたクラスとの距離の再計算のための準備
    num_h = element_class.eq(h).where.length
    num_k = element_class.eq(k).where.length
    dhk = dmin
    num_j = num_h+num_k

    #-- 融合された要素のクラス情報の更新
    class_num_new = element_class.max + 1
    [h,k].each{|cl_n|
      element_class[element_class.eq(cl_n)] = class_num_new
    }
    element_class_ary[true,istep+1] = element_class

    #-- 距離の計算
    class_list = element_class.to_a.uniq.sort
    nx = class_list.length
      deu_old = deu
      imax = (nx-1)*nx/2  # 組み合せの数 nP2
      deu = NArray.float(imax)
      oldcombclass = comb_class
      comb_class = Hash::new
      class_comb = Hash::new
    n=0
    (0..nx-2).each{|ni|
      (ni+1..nx-1).each{|nj|
        i = class_list[ni]
        j = class_list[nj]
        comb_class[[i,j]] = n
        class_comb[n] = [i,j]
        if [i,j].max < class_num_new
          #-- 融合以前から存在するクラスの保持
          deu[n] = deu_old[oldcombclass[[i,j]]]
        else
          #-- 新たに融合されたクラスとの距離の再計算
          num_i = element_class.eq(i).where.length
          dih = deu_old[oldcombclass[[i,h].sort]]
          dik = deu_old[oldcombclass[[i,k].sort]]
          deu[n] = ((num_i+num_h)*dih + (num_i+num_k)*dik - num_i*dhk) / (num_i+num_h+num_k)
        end
        n+=1
      }
    }

  nstep+=1
  }

=begin
  if renumber
    element_class_dum = element_class
    element_class_list = element_class.to_a.uniq
    i=0
    element_class_list.each{|el|
      element_class[element_class_dum.eq(el)] = i
      i+=1
    }
  end
=end
  if renumber
    element_class_dum = element_class
    element_class_list = element_class.to_a.uniq
    nclass = element_class_list.length
    cnum_list = NArray.int(nclass)
    i=0
    element_class_list.each{|el|
      cnum_list[i] = element_class_dum.eq(el).where.length
      i+=1
    }
    sort_index = cnum_list.sort_index
    i=0
    sort_index.each{|si|
      num = nclass-i-1
      el_dum = element_class_list[si]
      element_class[element_class_dum.eq(el_dum)] = num
      i+=1
    }
  end

  return [element_class,dmin_list]
  #return [element_class,dmin_list,element_class_ary,deu,class_comb]
end

#-----

if $0 == __FILE__

  data = NArray[[3.0,2.0,4.0,7.0,7.0,1.0,3.0]\
    ,[1.0,3.0,6.0,2.0,4.0,2.0,6.0]]
  nx,nt = data.shape

  #-- データの規格化
  data_stand = (data-data.mean)/data.stddev

  #-- クラスターアナリシス
  element_class, dmin_list \
    = cluster_analysis( data_stand, dim=0, nclass=3, renumber=true )
p element_class,dmin_list

  nx.times{|i|
    p [data[i,true],element_class[i]]
  }


  #--
begin
  require 'numru/ggraph'
  include NumRu

  require 'gphys_io_util'

  psname = $0.split('.')[0]

  ### === GGraph Sample Headers === eriko 29Jun2009 === ###
  DCL.swcset 'fname',psname

  iws = (ARGV[0]||(puts 'Workstation ID(I)?'; DCL.sgpwsn; gets)).to_i
  DCL.gropn iws

  DCLExt.sg_set_params 'lcntl'=>false#,'lfull'=>true
  DCLExt.uz_set_params 'indext1'=>3,'indext2'=>5\
    ,'indexl1'=>5,'indexl2'=>5,'inner'=>-1

  DCL.uzfact 0.6

  ### === GGraph Sample Headers === eriko 29Jun2009 === ###

  nx,nt = data.shape
  gphys1 = na2gphys(data[true,0],[NArray.sfloat(nx).indgen])
  gphys2 = na2gphys(data[true,1],[NArray.sfloat(nx).indgen])
  class_gphys = na2gphys(element_class,[NArray.sfloat(nx).indgen])

  GGraph.next_fig 'window'=>[0.0,10.0,0.0,10.0]
  #GGraph.scatter gphys1,gphys2,true#,'type'=>10
  GGraph.color_scatter gphys1,gphys2,class_gphys,true#,'type'=>10
  nx.times{|i|
    DCL.sgtxzu gphys1[i].val,gphys2[i].val,"#{element_class[i]}",0.02,0,0,5
  }

=begin
  GGraph.next_fig 'window'=>[0.0,3.0,0.0,1.5]
  GGraph.next_axes 'ytitle'=>'minimum distance between clusters'\
    ,'xside'=>'t','yunits'=>''

  nlength = dmin_list.length
  xax = NArray.sfloat(nlength).indgen
  gphys_dmin = na2gphys(dmin_list,[xax])
  GGraph.line gphys_dmin

  xlab=[]
  xax.each{|i| xlab << (nx-i).to_i.to_s }
  DCL.uxaxlb 'b',xax,xax,xlab,1
  DCL.uxsttl 'b','nums of clusters',0
=end

  DCL.grcls
rescue
end

end

