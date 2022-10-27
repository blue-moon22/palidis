/*
 * Install Interproscan
 */
process installInterproscan {

    input:


    output:
    path "${db_path}/${interproscan_db}"

    script:
    interproscan_link=params.interproscan_link
    tarball=params.interproscan_tarball
    interproscan_db=params.interproscan_db
    db_path=params.db_path

    """
    wget ${interproscan_link}
    wget ${interproscan_link}.md5
    md5sum -c ${tarball}.md5
    tar -pxvzf ${tarball}
    rm ${tarball}
    cd ${interproscan_db}

    # Edit interproscan.properties file
    sed -i 's/\${bin.directory}\\/prosite\\/pfscan/pfscan/' interproscan.properties
    sed -i 's/\${bin.directory}\\/prosite\\/pfsearch/pfsearch/' interproscan.properties
    sed -i 's/pfsearch_wrapper.py/\${bin.directory}\\/prosite\\/pfsearch_wrapper.py/' interproscan.properties

    python3 initial_setup.py
    cd ..

    mv ${interproscan_db} ${db_path}/
    """
}
