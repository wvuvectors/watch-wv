CREATE USER 'labuser'@'localhost' IDENTIFIED BY 'password123';
GRANT ALL PRIVILEGES ON patchr_db.* TO 'labuser'@'localhost';
FLUSH PRIVILEGES;