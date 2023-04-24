# strの中でtargetよりも前の部分文字列を返す
substr_str_before(str, target)=substr(str, 1, strstrt(str,target)-1)
# strの中でtargetよりも後の部分文字列を返す
substr_str_after(str, target)=\
     substr(str, strstrt(str,target)+strlen(target), strlen(str))
# strの中でindexより後の部分の文字列を返す
substr_after(str, index) = substr(str, index, strlen(str))
# strの中でindexより前の部分の文字列を返す
substr_before(str, index) = substr(str, 1, index)
# 再帰表現を利用して、strのindexより後の文字列のtargetをsubstに全て置換する
strsubst_sub(str, index, target, subst)=\
     strstrt(substr_after(str,index),target)==0 ? str :\
     strsubst_sub( substr_before(str, index-1) . substr_str_before(substr_after(str,index),target) . \
     subst . substr_str_after(substr_after(str,index),target),\
     index+strstrt(substr_after(str,index),target)+strlen(subst)-strlen(target),\
     target, subst)

# strの中のtargetをsubstに置換する
strsubst(str, target, subst)=strsubst_sub(str,1,target,subst)


#大きい方の値を返す関数を補助的に定義
max(x,y)= ( (x) > (y) ) ? (x) : (y)
# strのindexより後の文字列の中でtargetが最後に出てくる場所を返す
strstrlt_sub(str,index,target)=(strstrt(substr_after(str, index),target)==0 ? \
     max(index-strlen(target),0) : \
     strstrlt_sub(str, index-1+strstrt(substr_after(str, index),target)+strlen(target),target))
# strの中でtargetが最後に出てくる場所を返す
strstrlt(str,target)=strstrlt_sub(str,1,target)


# pathからディレクトリ名を返す。
dirname(path)=strstrt(path, "/")!= 0 ? substr(path, 1, strstrlt(path,"/")) : ""

# pathからディレクトリ名を消したものを返す。
remdirname(path)=substr(path,strstrlt(path,"/")+1,strlen(path))

# filenameからexentionを消したものを返す
remext(filename)=strstrt(filename, ".")==0 ? filename : \
                 substr(filename,1,strstrlt(filename, ".")-1)
