Do the same thing as too_elab, but have all the events compacted to within 2 years


ls *10yrs.db | xargs -I'{}' echo "glance_dir --db '{}'" > maf.sh
# ls *10yrs.db | xargs -I'{}' echo "scimaf_dir --db '{}'" >> maf.sh
ls *10yrs.db | xargs -I'{}' echo "scimaf_dir --limited --db '{}'" >> maf.sh
#ls *10yrs.db | xargs -I'{}' echo "ddf_dir --db '{}'" >> maf.sh
# ls *10yrs.db | xargs -I'{}' echo "metadata_dir --db '{}'" >> maf.sh

generate_ss 
cat ss_script.sh >> maf.sh

cat maf.sh | parallel -j 10
