from ntsparksearch.Common.Constants import Constants
import smtplib


class EmailSender(object):
    """Email model"""

    def __init__(self, receivers: list):
        self._receivers = receivers

    def _get_receivers(self):
        return self._receivers

    def mail_sender(self, msg):
        server = smtplib.SMTP(host=Constants.MAIL_HOST)
        server.ehlo()
        server.starttls()
        server.login(Constants.MAIL_USER, Constants.MAIL_PASS)
        server.sendmail(Constants.MAIL_SENDER, self._receivers, msg)

    def send_email_download_initialize(self, list_of_genes_to_download: list):
        myString = "\n".join(list_of_genes_to_download)

        if list_of_genes_to_download is not None:
            message_download = Constants.MSG_DOWNLOAD_INITIALIZE + myString
            try:
                self.mail_sender(message_download)
            except Exception as error:
                print("error" + repr(error))

    def send_email_download_finished(self, list_of_genes_to_download_finished: list):
        myString = "\n".join(list_of_genes_to_download_finished)

        if list_of_genes_to_download_finished is not None:
            message_download = Constants.MSG_DOWNLOAD_FINISHED + myString
            try:
                self.mail_sender(message_download)
            except Exception as error:
                print("error" + repr(error))