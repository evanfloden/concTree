while IFS=\= read var value; do
    vars+=($var)
    values+=($value)
done < nextflow.config

printf '%s\n' "${vars[@]}"
