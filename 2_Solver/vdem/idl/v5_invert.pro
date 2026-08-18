FUNCTION v5_invert,matrix0,data

  LUDC, matrix0, index, DOUBLE=double, INTERCHANGES=interchanges
  result = LUSOL(matrix0, index,data,/double)

return,result
END
